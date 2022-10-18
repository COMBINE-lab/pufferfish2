use kmers::naive_impl::CanonicalKmer;
use log::info;
use rand::{rngs::StdRng, Rng, SeedableRng};

use std::fmt::Debug;

use crate::index::SparseIndex;
use crate::index::{DecodeHit, Pufferfish, PufferfishBase};

pub mod refseq;

#[derive(Debug)]
pub struct SimpleSummary {
    pub n_kmers: u64,
    pub total_occs: u64,
    pub total_decoded: u64,
}

pub trait Benchmarking {
    fn decode_all_tps(&self, kmers: &[CanonicalKmer]) -> SimpleSummary;
}

impl<T> Benchmarking for T
where
    T: Pufferfish,
{
    fn decode_all_tps(&self, kmers: &[CanonicalKmer]) -> SimpleSummary {
        let mut decoded = 0;
        for km in kmers {
            let hits = self.get_ref_pos(km).unwrap(); // panics if not true positive;
            let iter = hits.iter();
            // let iter = DecodeHitIterator::new(&hits);
            let hits = Vec::from_iter(iter);
            decoded += hits.len();
        }

        SimpleSummary {
            n_kmers: kmers.len() as u64,
            total_occs: decoded as u64,
            total_decoded: decoded as u64,
        }
    }
}

//
// Random Kmer Generator Trait
//

pub trait RngKmer {
    fn rng_kmer_k(&self) -> u64;
    fn next_kmer_u64(&mut self) -> u64;
    fn next_kmer(&mut self) -> CanonicalKmer {
        let km = self.next_kmer_u64();
        CanonicalKmer::from_u64(km, self.rng_kmer_k() as u8)
    }
}

//
// Random Kmer generator
//

pub struct RandKmer<R> {
    rng: R,
    k: u64,
}

impl<R> RandKmer<R> {
    pub fn new(k: usize, rng: R) -> Self {
        Self { rng, k: k as u64 }
    }
}

impl<R: rand::RngCore> RngKmer for RandKmer<R> {
    fn rng_kmer_k(&self) -> u64 {
        self.k
    }

    fn next_kmer_u64(&mut self) -> u64 {
        self.rng.next_u64()
    }
}

//
// Random true positive kmer generator from index
//

pub struct IndexRng<'a, I: PufferfishBase, R> {
    rng: R,
    index: &'a I,
}

impl<'a, I: PufferfishBase, R> IndexRng<'a, I, R> {
    pub fn new(index: &'a I, rng: R) -> Self {
        Self { index, rng }
    }
}

impl<'a, I: PufferfishBase, R: rand::Rng> RngKmer for IndexRng<'a, I, R> {
    fn rng_kmer_k(&self) -> u64 {
        self.index.k() as u64
    }
    fn next_kmer_u64(&mut self) -> u64 {
        let n_refs = self.index.num_refs();
        let k = self.index.k();
        let ref_id = self.rng.gen_range(0..n_refs);
        let end_pos = self.index.ref_len(ref_id) - k;
        let start_pos = self.rng.gen_range(0..end_pos);

        self.index.get_refseq_kmer_u64(ref_id, start_pos)
    }
}

//
// Kmer bins
//

#[derive(Debug, Clone)]
pub struct KmerBinSlice<'a> {
    pub max_freq: Option<u64>,
    pub min_freq: Option<u64>,
    pub avg_freq: Option<f64>,

    k: u64,
    pub kmers: &'a [u64],
}

#[derive(Debug)]
pub struct KmerBin {
    pub max_freq: Option<u64>,
    pub min_freq: Option<u64>,
    pub avg_freq: Option<f64>,

    k: u64,
    pub kmers: Vec<u64>,
}

// Like the Index trait but instead of getting a T, you get a kmer
pub trait IndexKmer {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn index_kmer_k(&self) -> u64;
    fn get_kmer_u64(&self, i: usize) -> u64;
    fn get_kmer(&self, i: usize) -> CanonicalKmer {
        let km = self.get_kmer_u64(i);
        CanonicalKmer::from_u64(km, self.index_kmer_k() as u8)
    }
}

impl IndexKmer for KmerBinSlice<'_> {
    fn len(&self) -> usize {
        self.kmers.len()
    }

    fn index_kmer_k(&self) -> u64 {
        self.k
    }

    fn get_kmer_u64(&self, i: usize) -> u64 {
        self.kmers[i]
    }
}

impl IndexKmer for KmerBin {
    fn len(&self) -> usize {
        self.kmers.len()
    }

    fn index_kmer_k(&self) -> u64 {
        self.k
    }

    fn get_kmer_u64(&self, i: usize) -> u64 {
        self.kmers[i]
    }
}

pub struct KmerSampler<R: rand::Rng, K: IndexKmer> {
    rng: R,
    kmers: K,
}

// Genericeized KmerSampler, can sample a kmer given any Rng and IndexKmer
impl<R: rand::Rng, K: IndexKmer> KmerSampler<R, K> {
    pub fn new(rng: R, kmers: K) -> Self {
        Self { rng, kmers }
    }
}

impl<R: rand::Rng, K: IndexKmer> RngKmer for KmerSampler<R, K> {
    fn rng_kmer_k(&self) -> u64 {
        self.kmers.index_kmer_k()
    }
    fn next_kmer_u64(&mut self) -> u64 {
        let i = self.rng.gen_range(0..self.kmers.len());
        self.kmers.get_kmer_u64(i)
    }
}

//
// WIP: Binning unitigs and kmers by their frequency of occurence
//

#[derive(Debug)]
pub struct UnitigBins {
    pub max_freqs: Vec<u64>,
    pub avg_freqs: Vec<f64>,
    pub unitig_ids: Vec<Vec<u64>>,
}

#[derive(Debug)]
pub struct ExperimentKmers {
    pub kmers_per_bin: u64,
    pub max_freqs: Vec<u64>, // delimiters for bins
    pub avg_freqs: Vec<f64>,
    pub kmers_u64: Vec<Vec<u64>>,
    pub k: u64,
}

impl ExperimentKmers {
    pub fn n_bins(&self) -> usize {
        self.kmers_u64.len()
    }

    pub fn get_bin(&self, i: usize) -> KmerBinSlice {
        assert!(i < self.n_bins());
        let min_freq = if i == 0 { 0 } else { self.max_freqs[i - 1] };

        KmerBinSlice {
            min_freq: Some(min_freq),
            max_freq: Some(self.max_freqs[i]),
            avg_freq: Some(self.avg_freqs[i]),
            kmers: &self.kmers_u64[i],
            k: self.k,
        }
    }
}

// TODO impl this for all indices
pub trait BinKmers {
    fn bin_unitigs_by_freq(&self, n_bins: usize) -> UnitigBins;
    fn gen_kmers_binned_by_occ_freq(&self, n_bins: usize, kmers_per_bin: usize) -> ExperimentKmers;
}

impl BinKmers for SparseIndex {
    fn bin_unitigs_by_freq(&self, n_bins: usize) -> UnitigBins {
        let ctable = self.get_ctg_table();
        let offsets = ctable.get_offsets();

        // Get max contig freq
        let max_freq = {
            let mut last = 0;
            let mut max_freq = 0;
            for offset in offsets.iter().skip(1) {
                let n_occs = offset - last;
                max_freq = std::cmp::max(max_freq, n_occs);
                last = offset;
            }
            max_freq
        };

        info!("Max contig occ freq: {}", max_freq);

        let bin_width = (max_freq as f64).log2() / (n_bins as f64);
        let mut max_freqs: Vec<u64> = (0..n_bins)
            .into_iter()
            .map(|i| ((i + 1) as f64 * bin_width).exp2() as u64)
            .collect();

        max_freqs.push(max_freq);

        info!("freq delims {:?}", &max_freqs);

        let mut unitig_ids = vec![Vec::new(); n_bins + 1];
        let mut unitig_occs = vec![Vec::new(); n_bins + 1];

        let mut last = 0;
        for (ctg_id, offset) in offsets.iter().skip(1).enumerate() {
            let n_occs = offset - last;
            last = offset;

            let log_n_occs = (n_occs as f64).log2();

            let bin_i = (log_n_occs / bin_width) as u64;
            unitig_ids[bin_i as usize].push(ctg_id as u64);
            unitig_occs[bin_i as usize].push(n_occs as u64);
        }

        let avg_freqs: Vec<f64> = unitig_occs
            .into_iter()
            .map(|occs| {
                occs.iter()
                    .map(|&occ| occ as f64 / (occs.len() as f64))
                    .sum()
            })
            .collect();
        info!("max freqs {:?}", max_freqs);
        info!("avg freqs {:?}", avg_freqs);
        let lens: Vec<usize> = unitig_ids.iter().map(|b| b.len()).collect();
        info!("n_unitigs in bin {:?}", lens);

        // iterate over the contigs and bin them
        UnitigBins {
            avg_freqs,
            max_freqs,
            unitig_ids,
        }
    }

    fn gen_kmers_binned_by_occ_freq(&self, n_bins: usize, kmers_per_bin: usize) -> ExperimentKmers {
        let unitig_bins = self.bin_unitigs_by_freq(n_bins);

        let n_bins = unitig_bins.unitig_ids.len();

        let mut kmers_u64: Vec<Vec<u64>> = vec![Vec::with_capacity(kmers_per_bin); n_bins];

        let _min_freq = 0;

        for (bb, unitigs) in unitig_bins.unitig_ids.iter().enumerate() {
            // select a unitig
            let mut rng = StdRng::seed_from_u64(12091450981);
            let n_unitigs = unitigs.len();
            if n_unitigs > 0 {
                for _ in 0..kmers_per_bin {
                    let ii = rng.gen_range(0..n_unitigs);
                    let ctg_id = unitigs[ii] as usize;

                    let ctg_len = self.contig_len(ctg_id);
                    let offset = rng.gen_range(0..=(ctg_len - self.k()));
                    let km = self.get_contig_fw_kmer_u64(ctg_id, offset);
                    kmers_u64[bb].push(km);
                }
            }
        }

        ExperimentKmers {
            k: self.k() as u64,
            kmers_per_bin: kmers_per_bin as u64,
            max_freqs: unitig_bins.max_freqs,
            avg_freqs: unitig_bins.avg_freqs,
            kmers_u64,
        }
    }
}
