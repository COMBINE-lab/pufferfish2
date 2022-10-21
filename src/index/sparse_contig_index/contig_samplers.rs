use rand;
use rand::Rng;
use serde::{Deserialize, Serialize};
use simple_sds::bit_vector::BitVector;
use simple_sds::ops::Rank;
use simple_sds::raw_vector::{AccessRaw, RawVector};

use crate::index::{DecodeHit, Pufferfish};
use kmers::naive_impl::CanonicalKmer;

use clap::Subcommand;
use log::debug;

#[derive(Subcommand, Serialize, Deserialize, Debug, Clone, Copy, Eq, PartialEq)]
pub enum ContigSamplingStrategy {
    Minimal,
    // Random {
    //     #[clap(short, long)]
    //     p: f64,
    // },
    StridedBV {
        #[clap(short, long)]
        stride: usize,
    },

    StridedBVPop {
        #[clap(short, long)]
        stride: usize,

        #[clap(short, long)]
        thresh: usize,
    },

    Strided {
        #[clap(short, long)]
        stride: usize,
    },

    StridedOccThresh {
        // Sample with worst case 'stride' len, and sample occs > thresh occurences
        #[clap(short, long)]
        stride: usize,

        #[clap(short, long)]
        thresh: usize,
    },
}

impl ContigSamplingStrategy {
    pub fn get_is_sampled_bv<T: Pufferfish>(&self, pi: &T) -> BitVector {
        debug!("Getting is_sampled_bv with strategy {:?}", self);
        match *self {
            Self::Minimal => {
                let sampler = MinimalConditionsSampler::new(pi);
                sampler.get_is_sampled_bv()
            }
            Self::StridedBV { stride } => {
                let sampler = StridedBVSampler::new(pi, stride);
                sampler.get_is_sampled_bv()
            }
            Self::Strided { stride } => {
                let sampler = StridedSampler::new(pi, stride);
                sampler.get_is_sampled_bv()
            }
            Self::StridedOccThresh { stride, thresh } => {
                let sampler = StridedOccThreshSampler::new(pi, stride, thresh);
                sampler.get_is_sampled_bv()
            }
            Self::StridedBVPop { stride, thresh } => {
                let sampler = StridedBVPopSampler::new(pi, stride, thresh);
                sampler.get_is_sampled_bv()
            }
        }
    }
}

pub trait ContigSampler<'a, T> {
    fn get_is_sampled_bv(self) -> BitVector;
}

pub struct MinimalConditionsSampler<'a, T> {
    index: &'a T,
}

#[allow(dead_code)]
pub struct RandomContigSampler<'a, T> {
    index: &'a T,
    sampled_p: f64,
}

pub struct StridedBVPopSampler<'a, T> {
    // Samples every s-th index on the is_sampled_bv_vector and samples popular
    index: &'a T,
    stride: usize,
    thresh: usize,
}

pub struct StridedBVSampler<'a, T> {
    // Samples every s-th index on the is_sampled_bv_vector
    index: &'a T,
    stride: usize,
}

pub struct StridedSampler<'a, T> {
    // Samples every s-th index on reference tilings greedily
    index: &'a T,
    stride: usize,
}
pub struct StridedOccThreshSampler<'a, T> {
    // First samples contigs that occur > thresh times then
    // samples every s-th index on reference tilings greedily
    index: &'a T,
    stride: usize,
    thresh: usize,
}

impl<'a, T> MinimalConditionsSampler<'a, T> {
    pub fn new(index: &'a T) -> Self {
        Self { index }
    }
}

fn minimal_bv<T: Pufferfish>(index: &T) -> RawVector {
    let n_contigs = index.num_contigs();
    let mut bv = RawVector::with_len(n_contigs, false);

    for refseq in index.iter_refs() {
        let last_kmer_pos = refseq.ref_len() - index.k();
        let km_start = CanonicalKmer::from(refseq.get_kmer(0));
        let km_end = CanonicalKmer::from(refseq.get_kmer(last_kmer_pos));

        let start_ctg_id = index.get_ref_pos(&km_start).unwrap().contig_id();
        let end_ctg_id = index.get_ref_pos(&km_end).unwrap().contig_id();
        bv.set_bit(start_ctg_id, true);
        bv.set_bit(end_ctg_id, true);
    }
    bv
}

impl<'a, T> ContigSampler<'a, T> for MinimalConditionsSampler<'a, T>
where
    T: Pufferfish,
{
    fn get_is_sampled_bv(self) -> BitVector {
        let bv = minimal_bv(self.index);

        let mut is_sampled_bv = BitVector::from(bv);
        is_sampled_bv.enable_rank();
        is_sampled_bv
    }
}

#[allow(dead_code)]
impl<'a, T> RandomContigSampler<'a, T> {
    fn new(index: &'a T, p: f64) -> Self {
        assert!(p < 1., "Sampling prob must be < 1, {} given", p);
        Self {
            index,
            // nonsamp_ratio,
            sampled_p: p,
        }
    }
}

impl<'a, T> ContigSampler<'a, T> for RandomContigSampler<'a, T>
where
    T: Pufferfish,
{
    fn get_is_sampled_bv(self) -> BitVector {
        // TODO: could use rand::distributiosn::Bernoulli
        let mut bv = minimal_bv(self.index);
        for i in 0..bv.len() {
            let v = rand::thread_rng().gen_range(0.0..1.);
            if v < self.sampled_p {
                bv.set_bit(i, true);
            }
        }

        let mut is_sampled_bv = BitVector::from(bv);
        is_sampled_bv.enable_rank();
        is_sampled_bv
    }
}

impl<'a, T> StridedBVPopSampler<'a, T> {
    pub fn new(index: &'a T, stride: usize, thresh: usize) -> Self {
        Self {
            index,
            stride,
            thresh,
        }
    }
}

impl<'a, T> ContigSampler<'a, T> for StridedBVPopSampler<'a, T>
where
    T: Pufferfish,
{
    fn get_is_sampled_bv(self) -> BitVector {
        // TODO: could use rand::distributiosn::Bernoulli
        let mut bv = minimal_bv(self.index);
        let num_ctgs = self.index.num_contigs();

        for i in 0..num_ctgs {
            let n_occs = self.index.num_ctg_occs(i);
            let is_pop = n_occs >= self.thresh;
            let is_mod = i % self.stride == 0;
            if is_pop || is_mod {
                bv.set_bit(i, true);
            }
        }

        let mut is_sampled_bv = BitVector::from(bv);
        is_sampled_bv.enable_rank();
        is_sampled_bv
    }
}
impl<'a, T> StridedBVSampler<'a, T> {
    pub fn new(index: &'a T, stride: usize) -> Self {
        Self { index, stride }
    }
}

impl<'a, T> ContigSampler<'a, T> for StridedBVSampler<'a, T>
where
    T: Pufferfish,
{
    fn get_is_sampled_bv(self) -> BitVector {
        // TODO: could use rand::distributiosn::Bernoulli
        let mut bv = minimal_bv(self.index);
        for i in (0..bv.len()).step_by(self.stride) {
            bv.set_bit(i, true);
        }

        let mut is_sampled_bv = BitVector::from(bv);
        is_sampled_bv.enable_rank();
        is_sampled_bv
    }
}

impl<'a, T> StridedSampler<'a, T> {
    pub fn new(index: &'a T, stride: usize) -> Self {
        Self { index, stride }
    }
}

impl<'a, T> ContigSampler<'a, T> for StridedSampler<'a, T>
where
    T: Pufferfish,
{
    fn get_is_sampled_bv(self) -> BitVector {
        let mut bv = minimal_bv(self.index);

        // distance to last sampled contig
        let mut gap = 1;
        for refseq in self.index.iter_refs() {
            debug!("ref id: {}", refseq.ref_id);
            for ctg in refseq.iter_contigs() {
                let is_sampled = bv.bit(ctg.id);
                if is_sampled {
                    // reset the gap
                    gap = 1;

                    continue;
                } else if gap >= self.stride {
                    // reset the gap and sample
                    gap = 1;
                    bv.set_bit(ctg.id, true);
                } else {
                    gap += 1;
                }
            }
        }

        let mut is_sampled_bv = BitVector::from(bv);
        is_sampled_bv.enable_rank();
        is_sampled_bv
    }
}

impl<'a, T> StridedOccThreshSampler<'a, T> {
    pub fn new(index: &'a T, stride: usize, thresh: usize) -> Self {
        Self {
            index,
            stride,
            thresh,
        }
    }
}

impl<'a, T> ContigSampler<'a, T> for StridedOccThreshSampler<'a, T>
where
    T: Pufferfish,
{
    fn get_is_sampled_bv(self) -> BitVector {
        // TODO: could use rand::distributiosn::Bernoulli
        let mut bv = minimal_bv(self.index);

        // distance to last sampled contig
        let mut gap = 0;
        for refseq in self.index.iter_refs() {
            for ctg in refseq.iter_contigs() {
                let is_sampled = bv.bit(ctg.id);
                if is_sampled {
                    // reset the gap
                    gap = 0;
                    continue;
                } else if (ctg.num_occs > self.thresh) || (gap >= self.stride) {
                    // reset the gap and sample
                    gap = 0;
                    bv.set_bit(ctg.id, true);
                } else {
                    gap += 1;
                }
            }
        }

        let mut is_sampled_bv = BitVector::from(bv);
        is_sampled_bv.enable_rank();
        is_sampled_bv
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::index::sparse_contig_index::WMSparseSparseIndex;
    use crate::index::validate::*;
    use crate::index::{DenseIndex, SparseIndex};
    use crate::test_utils::*;

    #[test]
    fn strided_bv() {
        //TODO move this to integration tests,
        // and impl Pufferfish for SparseContigdenseindex
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        let strategy = ContigSamplingStrategy::StridedBV { stride: 3 };
        let pi = WMSparseSparseIndex::from_sparse_index_w_strategy(pi, strategy, 0, 0);
        assert_is_valid(&pi);

        assert_eq!(pi.num_sampled_contigs(), 18);
    }

    #[test]
    fn strided() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        // Two references with the following contig tilings:
        // minimal sampler samples first and last contig
        // 0-1-2-3-4
        // 0-5-2-6-4

        let sampler = StridedSampler::new(&pi, 1);
        let bv = sampler.get_is_sampled_bv();
        let bv_true = RawVector::with_len(7, true);
        assert_eq!(RawVector::from(bv), bv_true);

        let sampler = StridedSampler::new(&pi, 2);
        let bv = sampler.get_is_sampled_bv();
        let mut bv_true = RawVector::with_len(7, false);
        bv_true.set_bit(0, true);
        bv_true.set_bit(2, true);
        bv_true.set_bit(4, true);
        assert_eq!(RawVector::from(bv), bv_true);

        let sampler = StridedSampler::new(&pi, 3);
        let bv = sampler.get_is_sampled_bv();
        let mut bv_true = RawVector::with_len(7, false);
        bv_true.set_bit(0, true);
        bv_true.set_bit(3, true);
        bv_true.set_bit(6, true);
        bv_true.set_bit(4, true);
        assert_eq!(RawVector::from(bv), bv_true);

        let sampler = StridedSampler::new(&pi, 4);
        let bv = sampler.get_is_sampled_bv();
        let mut bv_true = RawVector::with_len(7, false);
        bv_true.set_bit(0, true);
        bv_true.set_bit(4, true);
        assert_eq!(RawVector::from(bv), bv_true);

        let sampler = StridedSampler::new(&pi, 1000000);
        let bv = sampler.get_is_sampled_bv();
        let mut bv_true = RawVector::with_len(7, false);
        bv_true.set_bit(0, true);
        bv_true.set_bit(4, true);
        assert_eq!(RawVector::from(bv), bv_true);
    }
}
