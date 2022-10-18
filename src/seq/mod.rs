use serde::{Deserialize, Serialize};
use simple_sds::int_vector::IntVector;
use simple_sds::ops::Vector;
use simple_sds::raw_vector::{AccessRaw, RawVector};

use crate::serde_ext;
use kmers::naive_impl::Kmer;

#[derive(Clone, Debug, Serialize, Deserialize, Eq, PartialEq)]
pub struct SeqVector {
    // k: u8,
    #[serde(with = "serde_ext")]
    data: RawVector,
}

#[derive(Clone)]
pub struct SeqVectorSlice<'a> {
    len: usize,
    start_pos: usize,
    slice: &'a SeqVector, // TODO implement RawVectorSlice? instead of pointing
}

impl SeqVectorSlice<'_> {
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get_kmer(&self, pos: usize, k: u8) -> Kmer {
        let km = self.get_kmer_u64(pos, k);
        Kmer::from_u64(km, k)
    }

    pub fn get_kmer_u64(&self, pos: usize, k: u8) -> u64 {
        assert!(pos < self.len());
        let pos = pos + self.start_pos;
        self.slice.get_kmer_u64(pos, k)
    }

    pub fn get_base(&self, pos: usize) -> u64 {
        self.get_kmer_u64(pos, 1)
    }

    pub fn slice(&self, start: usize, end: usize) -> Self {
        assert!(end <= self.len());
        Self {
            len: end - start,
            start_pos: start,
            slice: self.slice,
        }
    }

    pub fn iter_kmers(&self, k: u8) -> SeqVecKmerIterator {
        SeqVecKmerIterator {
            k,
            len: self.len - (k as usize) + 1,
            pos: 0,
            seq: self.clone(),
        }
    }
}

impl SeqVector {
    pub fn len(&self) -> usize {
        self.data.len() / 2
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn get_kmer(&self, pos: usize, k: u8) -> Kmer {
        Kmer::from_u64(self.get_kmer_u64(pos, k), k)
    }

    pub fn get_kmer_u64(&self, pos: usize, k: u8) -> u64 {
        assert!(pos < self.len());
        unsafe { self.data.int(pos * 2, k as usize * 2) }
    }

    pub fn get_base(&self, pos: usize) -> u64 {
        self.get_kmer_u64(pos, 1)
    }

    pub fn as_slice(&self) -> SeqVectorSlice<'_> {
        SeqVectorSlice {
            start_pos: 0,
            len: self.len(),
            slice: self,
        }
    }

    pub fn slice(&self, start: usize, end: usize) -> SeqVectorSlice {
        self.as_slice().slice(start, end)
    }

    pub fn iter_kmers(&self, k: u8) -> SeqVecKmerIterator {
        SeqVecKmerIterator {
            k,
            len: self.len() - (k as usize) + 1,
            pos: 0,
            seq: self.as_slice(),
        }
    }
}

impl std::fmt::Display for SeqVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SeqVector[ {} ]", String::from(self))
    }
}

impl From<&SeqVector> for String {
    fn from(data: &SeqVector) -> Self {
        let mut str = String::new();
        let bases = vec!['A', 'C', 'G', 'T'];
        for i in 0..data.len() {
            let base = data.get_base(i);
            let base = bases[base as usize];
            str.push(base);
        }
        str
    }
}
impl From<SeqVector> for String {
    fn from(data: SeqVector) -> Self {
        Self::from(&data)
    }
}

impl From<&String> for SeqVector {
    fn from(data: &String) -> Self {
        assert!(data.is_ascii());
        let bytes = data.as_bytes();
        Self::from(bytes)
    }
}

impl From<String> for SeqVector {
    fn from(data: String) -> Self {
        Self::from(&data)
    }
}

impl From<&[u8]> for SeqVector {
    fn from(data: &[u8]) -> Self {
        let len = data.len() * 2;
        let chunks = data.chunks(32);
        let mut words = Vec::with_capacity(chunks.len());
        for chunk in chunks {
            let word = Kmer::from(chunk);
            words.push(word.into_u64());
        }
        let rv = RawVector::from_parts(len, words);
        Self { data: rv }
    }
}

impl From<RawVector> for SeqVector {
    fn from(data: RawVector) -> Self {
        assert_eq!(data.len() % 2, 0);
        Self { data }
    }
}

impl From<IntVector> for SeqVector {
    fn from(data: IntVector) -> Self {
        assert_eq!(data.width(), 2);
        Self {
            data: RawVector::from(data),
        }
    }
}

pub struct SeqVecKmerIterator<'a> {
    k: u8,
    len: usize,
    pos: usize,
    seq: SeqVectorSlice<'a>,
}

impl<'a> SeqVecKmerIterator<'a> {
    pub fn new(slice: SeqVectorSlice<'a>, k: u8) -> Self {
        Self {
            k,
            len: slice.len() - (k as usize) + 1,
            pos: 0,
            seq: slice,
        }
    }
}

impl SeqVecKmerIterator<'_> {
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Iterator for SeqVecKmerIterator<'_> {
    type Item = Kmer;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos < self.len() {
            let km = self.seq.get_kmer(self.pos, self.k);
            self.pos += 1;
            Some(km)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn seq_slice_test() {
        let bytes = vec![1u64, 2, 3];
        let iv = IntVector::from(bytes);
        let rv = RawVector::from(iv);
        let sv = SeqVector::from(rv);

        let slice = sv.as_slice();

        assert_eq!(slice.len(), 32 * 3);
        assert_eq!(slice.get_kmer_u64(0, 32), 1);

        let slice = sv.slice(1, 96);
        assert_eq!(slice.get_kmer_u64(0, 32), sv.get_kmer_u64(1, 32));

        let slice = sv.slice(75, 96);
        assert_eq!(slice.get_kmer_u64(0, 7), sv.get_kmer_u64(75, 7));
    }
}
