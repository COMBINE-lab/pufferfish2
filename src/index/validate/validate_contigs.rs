use kmers::naive_impl::CanonicalKmer;
use serde::{Deserialize, Serialize};

use crate::index::{CachedQuery, DecodeHit, MappedRefPos, Pufferfish, QueryCache};

use log::debug;

pub trait GenContigMRPs {
    fn gen_validation_kmers(&self, debug: usize) -> ContigKmerMRPs;
}

#[derive(Serialize, Deserialize)]
pub struct ContigKmerMRPs {
    seq_hash: String,
    k: usize,
    first_kmers: Vec<u64>,
    last_kmers: Vec<u64>,
    mid_kmers: Vec<u64>,
    contig_lens: Vec<usize>,

    first_mrps: Vec<MappedRefPos>,
    offsets: Vec<usize>,
}

impl<T> GenContigMRPs for T
where
    T: Pufferfish + for<'a> CachedQuery<'a, QueryCache>,
{
    fn gen_validation_kmers(&self, debug: usize) -> ContigKmerMRPs {
        let n_contigs = self.num_contigs();

        let mut first_kmers = Vec::with_capacity(n_contigs);
        let mut last_kmers = Vec::with_capacity(n_contigs);
        let mut mid_kmers = Vec::with_capacity(n_contigs);
        let mut offsets = Vec::with_capacity(n_contigs);
        let mut contig_lens = Vec::with_capacity(n_contigs);

        let mut first_mrps = Vec::new();
        let k = self.k() as usize;

        let mut offset = 0;

        // iterate over the contigs
        debug!(
            "Generating MRPs for prefixes, suffixes, and middle kmers for {} contigs",
            n_contigs
        );

        let debug = if debug > 0 { n_contigs / debug } else { 0 };
        // let mut cache = QueryCache::new(); // TODO: maybe use cache instead to be fast...
        for ctg_id in 0..n_contigs {
            if debug > 0 && ctg_id % debug == 0 {
                debug!(
                    "\t({:.1}%) {} / {}",
                    (ctg_id as f64) / (n_contigs as f64) * 100.0,
                    ctg_id,
                    n_contigs
                );
            }
            let ctg_len = self.contig_len(ctg_id);
            let prefix = self.get_contig_fw_kmer_u64(ctg_id, 0);
            let suffix: u64 = self.get_contig_fw_kmer_u64(ctg_id, ctg_len - k);
            let mid_kmer: u64 = self.get_contig_fw_kmer_u64(ctg_id, (ctg_len - k) / 2);

            contig_lens.push(ctg_len);
            first_kmers.push(prefix);
            last_kmers.push(suffix);
            mid_kmers.push(mid_kmer);

            let km = CanonicalKmer::from_u64(prefix, k as u8);

            // let mut mrps = self.get_ref_pos_with_cache(&km, &mut cache).unwrap().as_vec();
            let mut mrps = self.get_ref_pos(&km).unwrap().as_vec();

            let n_occs = mrps.len();
            offsets.push(offset);

            offset += n_occs;
            first_mrps.append(&mut mrps); // Last thing to do before mutating appended mrps just to be safe
        }
        offsets.push(offset);

        ContigKmerMRPs {
            seq_hash: self.get_base_index().seq_hash.clone(),
            k,
            first_kmers,
            last_kmers,
            mid_kmers,
            contig_lens,
            first_mrps,
            offsets,
        }
    }
}

impl ContigKmerMRPs {
    pub fn n_contigs(&self) -> usize {
        self.contig_lens.len()
    }

    pub fn n_occs(&self, ctg_id: usize) -> usize {
        self.offsets[ctg_id + 1] - self.offsets[ctg_id]
    }

    pub fn validate_index<T: Pufferfish>(&self, index: &T) {
        let debug = 20;
        assert_eq!(self.k, index.k());
        assert_eq!(self.seq_hash, index.get_base_index().seq_hash);

        let debug = self.n_contigs() / debug;
        let n_contigs = self.n_contigs();

        let now = std::time::Instant::now();
        for i in 0..n_contigs {
            if debug > 0 && i % debug == 0 {
                log::debug!(
                    "({:.1}%) {} / {}",
                    (i as f64) / (self.n_contigs() as f64) * 100.0,
                    i,
                    n_contigs
                )
            }
            let prefix = CanonicalKmer::from_u64(self.first_kmers[i], self.k as u8);
            let suffix = CanonicalKmer::from_u64(self.last_kmers[i], self.k as u8);
            let mid_kmer = CanonicalKmer::from_u64(self.mid_kmers[i], self.k as u8);
            let ctg_len = self.contig_lens[i];
            let n_occs = self.n_occs(i);
            // Check contig prefix MRPs

            let mrps = index.get_ref_pos(&prefix).unwrap();
            assert_eq!(mrps.contig_id(), i);
            assert_eq!(mrps.contig_len(), ctg_len);
            assert_eq!(mrps.len(), n_occs);
            let mrps = mrps.as_vec();
            self.check_ctg_mrp(i, &mrps, ctg_len, 0);

            let mrps = index.get_ref_pos(&suffix).unwrap();
            assert_eq!(mrps.contig_id(), i);
            assert_eq!(mrps.contig_len(), ctg_len);
            assert_eq!(mrps.len(), n_occs);
            let mrps = mrps.as_vec();
            self.check_ctg_mrp(i, &mrps, ctg_len, ctg_len - self.k);

            let mrps = index.get_ref_pos(&mid_kmer).unwrap();
            assert_eq!(mrps.contig_id(), i);
            assert_eq!(mrps.contig_len(), ctg_len);
            assert_eq!(mrps.len(), n_occs);
            let mrps = mrps.as_vec();
            self.check_ctg_mrp(i, &mrps, ctg_len, (ctg_len - self.k) / 2);
        }

        let elapsed = now.elapsed();
        let ns_per_kmer = elapsed.as_nanos() / (n_contigs as u128) / 3;
        let throughput = (n_contigs as f64) / elapsed.as_secs_f64();
        debug!("Validated kmers on {} contigs", n_contigs);
        debug!("Query time per kmer: {} ns", ns_per_kmer);
        debug!("Throughput: {:.2} kmers/sec", throughput);
    }

    fn check_ctg_mrp(
        &self,
        ctg_id: usize,
        mrps_to_check: &[MappedRefPos],
        ctg_len: usize,
        offset: usize,
    ) {
        let s = self.offsets[ctg_id];
        let e = self.offsets[ctg_id + 1];
        let valid_mrps = &self.first_mrps[s..e];

        for (check, valid) in mrps_to_check.iter().zip(valid_mrps) {
            let valid_pos = if valid.is_forward() {
                valid.pos + offset
            } else {
                valid.pos - offset
            };
            let is_ok = (check.orientation == valid.orientation)
                && (check.pos == valid_pos)
                && (check.ref_id == valid.ref_id);
            // assert_eq!(check.pos, valid_pos);
            if !is_ok {
                let valid = MappedRefPos {
                    pos: valid_pos,
                    ..*valid
                };
                panic!("Unitig {} w/ len {} validation error at offset {}.\n\tValid: {:?}\n\tFound: {:?}", ctg_id, ctg_len, offset, valid, check)
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::index::{DenseIndex, SparseIndex};
    use crate::test_utils::*;

    #[test]
    fn test() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        let ctg_kmers = pi.gen_validation_kmers(0);
        ctg_kmers.validate_index(&pi);

        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let ctg_kmers = pi.gen_validation_kmers(0);
        ctg_kmers.validate_index(&pi);

        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let spi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        ctg_kmers.validate_index(&spi);
    }
}
