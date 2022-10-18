use std::io::Result;
use std::path::Path;

use log::debug;
use serde::{Deserialize, Serialize};

// TODO... wtf are these paths.
use crate::index::BaseIndex;
use crate::seq::SeqVector;
use crate::{cpp::FromCereal, index::compact::FromCompact};

use super::IndexKmer;

#[derive(PartialEq, Debug, Clone, Serialize, Deserialize)]
pub struct RefSeq {
    ref_seq: SeqVector,
    ref_lens: Vec<u32>,
    ref_accum_lens: Vec<u64>,
}

impl RefSeq {
    pub fn deserialize_from_cpp<P: AsRef<Path>>(
        ref_seq: P,
        ref_lens: P,
        ref_accum_lens: P,
    ) -> Result<Self> {
        debug!("Loading refseq from {}", ref_seq.as_ref().to_string_lossy());
        let ref_seq = SeqVector::from_compact_serialized(ref_seq).unwrap();
        debug!(
            "Loading reflens from {}",
            ref_lens.as_ref().to_string_lossy()
        );
        //TODO assert that the len of ref_seq is the sum of lens of reference sequences
        let ref_lens: Vec<u32> = Vec::load_from_cereal_archive(ref_lens)?;

        debug!(
            "Loading ref_accum_lens from {}",
            ref_accum_lens.as_ref().to_string_lossy()
        );

        let ref_accum_lens: Vec<u64> = {
            // TODO: note, we are appending one zero to the front so accessing the offset is easy
            let mut v = Vec::load_from_cereal_archive(ref_accum_lens)?;
            v.insert(0, 0);
            v
        };

        Ok(Self {
            ref_seq,
            ref_lens,
            ref_accum_lens,
        })
    }

    pub fn num_refs(&self) -> usize {
        self.ref_lens.len()
    }

    pub fn tot_len(&self) -> usize {
        self.ref_accum_lens[self.num_refs()] as usize
    }

    pub fn num_kmers(&self, k: usize) -> usize {
        for &l in &self.ref_lens {
            assert!(l >= (k as u32));
        }

        self.tot_len() - (self.num_refs() * (k - 1))
    }
}

pub struct RefSeqKmers {
    ref_seq: RefSeq,
    partition_points: Vec<usize>,
    k: usize,
}

impl IndexKmer for RefSeqKmers {
    fn len(&self) -> usize {
        let last = self.partition_points.len() - 1;
        self.partition_points[last]
    }

    fn get_kmer_u64(&self, i: usize) -> u64 {
        let ref_id = self.partition_points.partition_point(|&p| i >= p);
        let offset = i + (ref_id * (self.k - 1));
        self.ref_seq.ref_seq.get_kmer_u64(offset, self.k as u8)
    }

    fn index_kmer_k(&self) -> u64 {
        self.k as u64
    }
}

impl RefSeqKmers {
    pub fn new(ref_seq: RefSeq, k: usize) -> Self {
        let partition_points = ref_seq.ref_accum_lens[1..]
            .iter()
            .enumerate()
            .map(|(i, &l)| (l as usize) - ((i + 1) * (k - 1)))
            .collect();
        Self {
            ref_seq,
            partition_points,
            k,
        }
    }
}

impl From<BaseIndex> for RefSeq {
    fn from(index: BaseIndex) -> Self {
        RefSeq {
            ref_seq: index.ref_seq.expect("No refseq found"),
            ref_lens: index.ref_lens,
            ref_accum_lens: index.ref_accum_lens,
        }
    }
}

impl From<&BaseIndex> for RefSeq {
    fn from(index: &BaseIndex) -> Self {
        RefSeq {
            ref_seq: index.ref_seq.clone().expect("No refseq found"),
            ref_lens: index.ref_lens.clone(),
            ref_accum_lens: index.ref_accum_lens.clone(),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::index::prelude::*;
    use crate::test_utils::*;
    #[test]
    fn ser_de_refseq() {
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let pi = pi.get_base_index();

        let ref_seq = RefSeq::from(pi);

        let ref_seq_from_file = RefSeq::deserialize_from_cpp(
            p.join("refseq.bin"),
            p.join("reflengths.bin"),
            p.join("refAccumLengths.bin"),
        )
        .unwrap();

        assert_eq!(ref_seq, ref_seq_from_file);
    }

    #[test]
    fn refseq_kmers() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let ref_seq = RefSeq::deserialize_from_cpp(
            p.join("refseq.bin"),
            p.join("reflengths.bin"),
            p.join("refAccumLengths.bin"),
        )
        .unwrap();

        assert_eq!(ref_seq.num_refs(), 2);
        assert_eq!(ref_seq.tot_len(), 19 + 21);
        assert_eq!(ref_seq.num_kmers(19), 1 + 3);
        assert_eq!(ref_seq.num_kmers(3), 17 + 19);
    }

    #[test]
    fn refseq_kmers_sampler_small() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let ref_seq = RefSeq::deserialize_from_cpp(
            p.join("refseq.bin"),
            p.join("reflengths.bin"),
            p.join("refAccumLengths.bin"),
        )
        .unwrap();

        let sampler = RefSeqKmers::new(ref_seq.clone(), 19);
        assert_eq!(sampler.partition_points, vec![1, ref_seq.num_kmers(19)]);

        let sampler = RefSeqKmers::new(ref_seq.clone(), 3);
        assert_eq!(sampler.partition_points, vec![17, ref_seq.num_kmers(3)]);

        // 0th kmer in ref 0: AGT
        let km = sampler.get_kmer_u64(0);
        println!("{:b}", km);
        let kw = 0b11_10_00;
        assert_eq!(km, kw);

        // 3rd kmer in ref 0: GAT
        let km = sampler.get_kmer_u64(3);
        let kw = 0b11_00_10;
        assert_eq!(km, kw);

        // last kmer in ref 0: GTA
        let km = sampler.get_kmer_u64(16);
        let kw = 0b00_11_10;
        assert_eq!(km, kw);

        // 0th kmer in ref 1: AGT
        let km = sampler.get_kmer_u64(17);
        let kw = 0b11_10_00;
        assert_eq!(km, kw);
    }
}
