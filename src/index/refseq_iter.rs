// use super::base_index::BaseIndex;
use kmers::naive_impl::{CanonicalKmer, Kmer};

use super::{ContigOrientation, DecodeHit, MappedOrientation, Pufferfish, PufferfishBase};
use crate::seq::*;

// Wrapper for a "reference" in an base index
pub struct BaseIndexReference<'a, T> {
    pub ref_id: usize,
    index: &'a T,
}

pub struct RefSeqContigIterator<'a, T> {
    // start position of current contig
    pos: usize,
    ref_id: usize,
    end_pos: usize,
    index: &'a T,
}

// Iterator over references
pub struct RefIterator<'a, T> {
    ref_id: usize,
    num_refs: usize,
    index: &'a T,
}

impl<T: PufferfishBase> BaseIndexReference<'_, T> {
    pub fn ref_len(&self) -> usize {
        self.index.ref_len(self.ref_id)
    }

    pub fn get_kmer(&self, pos: usize) -> Kmer {
        self.index.get_refseq_kmer(self.ref_id, pos)
    }
}

impl<'a, T: PufferfishBase> BaseIndexReference<'a, T> {
    pub fn iter_kmers(&self) -> SeqVecKmerIterator<'a> {
        // Iterate over the kmers of a reference
        let k = self.index.k() as u8;
        let seq = self.index.get_refseq(self.ref_id);
        SeqVecKmerIterator::new(seq, k)
    }

    pub fn iter_contigs(&self) -> RefSeqContigIterator<'a, T> {
        // Iterate over the contigs that tile a reference
        RefSeqContigIterator::new(self)
    }
}

impl<'a, T: PufferfishBase> RefIterator<'a, T> {
    pub fn new(index: &'a T) -> Self {
        Self {
            ref_id: 0,
            num_refs: index.num_refs(),
            index,
        }
    }
}

impl<'a, T: PufferfishBase> RefSeqContigIterator<'a, T> {
    pub fn new(ref_iter: &BaseIndexReference<'a, T>) -> Self {
        let ref_id = ref_iter.ref_id;
        let index = ref_iter.index;
        let end_pos = index.ref_len(ref_id) - (index.k() as usize) + 1;
        Self {
            pos: 0,
            end_pos,
            ref_id,
            index,
        }
    }
}

impl<'a, T: PufferfishBase> Iterator for RefIterator<'a, T> {
    type Item = BaseIndexReference<'a, T>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.ref_id < self.num_refs {
            let item = Self::Item {
                ref_id: self.ref_id,
                index: self.index,
            };
            self.ref_id += 1;
            Some(item)
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct RefSeqContigOcc {
    pub o: ContigOrientation,
    pub pos: usize,
    pub len: usize,
    pub id: usize,
    // Total number of occurences in the reference
    pub num_occs: usize,
}

impl<'a, T: Pufferfish> Iterator for RefSeqContigIterator<'a, T> {
    type Item = RefSeqContigOcc;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos < self.end_pos {
            let km = self.index.get_refseq_kmer(self.ref_id, self.pos);
            let km = CanonicalKmer::from(km);
            let hits = self.index.get_ref_pos(&km).unwrap();
            let o = match hits.mapped_orientation() {
                MappedOrientation::Forward => ContigOrientation::Forward,
                MappedOrientation::Backward => ContigOrientation::Backward,
            };

            let ctg = RefSeqContigOcc {
                o,
                pos: self.pos,
                len: hits.contig_len(),
                id: hits.contig_id(),
                num_occs: hits.len(),
            };

            self.pos += ctg.len - self.index.k() + 1;

            Some(ctg)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::index::DenseIndex;
    use crate::test_utils::*;

    #[test]
    fn test_contig_iter() {
        // Check unitig tiling of small example
        // See: /test_data/tiny-multi-refs
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        let refs = Vec::from_iter(pi.iter_refs());

        let ctgs0 = Vec::from_iter(refs[0].iter_contigs());
        let ctgs1 = Vec::from_iter(refs[1].iter_contigs());

        let lens0: Vec<usize> = ctgs0.iter().map(|x| x.len).collect();
        let occs0: Vec<usize> = ctgs0.iter().map(|x| x.num_occs).collect();
        let ids0: Vec<usize> = ctgs0.iter().map(|x| x.id).collect();

        let lens1: Vec<usize> = ctgs1.iter().map(|x| x.len).collect();
        let occs1: Vec<usize> = ctgs1.iter().map(|x| x.num_occs).collect();
        let ids1: Vec<usize> = ctgs1.iter().map(|x| x.id).collect();

        assert_eq!(ctgs0.len(), 5);
        assert_eq!(lens0, vec![5, 8, 9, 8, 5]);
        assert_eq!(occs0, vec![2, 1, 2, 1, 2]);
        assert_eq!(ids0, vec![0, 1, 2, 3, 4]);

        assert_eq!(ctgs1.len(), 5);
        assert_eq!(lens1, vec![5, 9, 9, 9, 5]);
        assert_eq!(occs1, vec![2, 1, 2, 1, 2]);
        assert_eq!(ids1, vec![0, 5, 2, 6, 4]);
    }

    #[test]
    fn test_refseq_kmer_iter() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        let ref1 = pi.get_refseq(1);

        assert_eq!(ref1.len(), 21);

        let k = pi.k() as u8;

        let kmers = vec![
            "agtga", "gtgac", "tgact", "gactg", "actga", "ctgat", "tgata", "gatag", "atagt",
            "tagta", "agtag", "gtagc", "tagca", "agcag", "gcagg", "caggt", "aggta",
        ];

        let kmers_ = Vec::from_iter(ref1.iter_kmers(k));
        let kmers_: Vec<String> = kmers_.iter().map(|km| format!("{}", km)).collect();
        assert_eq!(kmers, kmers_);
    }
}
