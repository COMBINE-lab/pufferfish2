use super::{
    CachedQuery, DecodeHit, MappedOrientation, MappedRefPos, PuffQuery, Pufferfish, QueryCache,
};
use kmers::naive_impl::CanonicalKmer;
use log::debug;

pub mod validate_contigs;

pub trait Validate {
    fn validate(&self) -> bool;
}

pub trait ValidateCached<'a, C> {
    fn validate_w_cache(&'a self, cache: &mut C) -> bool;
}

impl<T> Validate for T
where
    T: Pufferfish,
{
    fn validate(&self) -> bool {
        if !self.has_ref_seq() {
            // TODO message to tell user to get set pufferfish flag for ref seq
            panic!("Cannot validate index without stored reference sequence");
        }

        let n_refs = self.num_refs();

        for (ref_id, refseq) in self.iter_refs().enumerate() {
            debug!(
                "Validating reference {} of {}, with length {}",
                ref_id,
                n_refs,
                refseq.ref_len()
            );

            for (pos, kmer) in refseq.iter_kmers().enumerate() {
                let canon_kmer = CanonicalKmer::from(kmer);
                let mut is_ok = false;

                let hits = self.get_ref_pos(&canon_kmer).unwrap();

                let mrp = MappedRefPos {
                    ref_id,
                    pos,
                    orientation: MappedOrientation::Forward,
                };

                for mapped_pos in hits.as_vec() {
                    // TODO: could use hits.iter() to early exit
                    // but for some impls of DecodeHit, as_vec is super fast
                    let found = mapped_pos == mrp;
                    is_ok |= found;
                }

                if !is_ok {
                    panic!(
                        "WARN: validation failed @ ref: {}, kmer: {} {}",
                        ref_id, pos, canon_kmer
                    );
                }
            }
        }
        true
    }
}

impl<'a, T, C> ValidateCached<'a, C> for T
where
    T: Pufferfish + CachedQuery<'a, C>,
{
    fn validate_w_cache(&'a self, cache: &mut C) -> bool {
        if !self.has_ref_seq() {
            // TODO message to tell user to get set pufferfish flag for ref seq
            panic!("Cannot validate index without stored reference sequence");
        }

        // let mut cache = QueryCache::new();

        let n_refs = self.num_refs();

        for (ref_id, refseq) in self.iter_refs().enumerate() {
            debug!(
                "Validating reference {} of {}, with length {}",
                ref_id,
                n_refs,
                refseq.ref_len()
            );

            for (pos, kmer) in refseq.iter_kmers().enumerate() {
                let canon_kmer = CanonicalKmer::from(kmer);
                let mut is_ok = false;

                // Like validate, but also checks that cached and uncached queries return the same result
                // TODO check contig IDs and offsets into contigs are also correct.
                let hits_not_cached = self.get_ref_pos(&canon_kmer).unwrap();
                let hits = self.get_ref_pos_with_cache(&canon_kmer, cache).unwrap();

                let mrps_not_cached = hits_not_cached.as_vec();
                let mrps_cached = hits.as_vec();
                assert_eq!(mrps_cached, mrps_not_cached);

                let mrp = MappedRefPos {
                    ref_id,
                    pos,
                    orientation: MappedOrientation::Forward,
                };

                for mapped_pos in &mrps_cached {
                    // TODO: could use hits.iter() to early exit
                    // but for some impls of DecodeHit, as_vec is super fast
                    let found = mapped_pos == &mrp;
                    is_ok |= found;
                }

                if !is_ok {
                    panic!(
                        "WARN: validation failed @ ref: {}, kmer: {} {}",
                        ref_id, pos, canon_kmer
                    );
                }
            }
        }
        true
    }
}

pub fn assert_is_valid<T>(v: &T)
where
    T: Validate,
{
    let is_valid = v.validate();
    assert!(is_valid)
}

pub fn assert_is_valid_cached<'a, T>(v: &'a T)
where
    T: ValidateCached<'a, QueryCache>,
{
    let mut cache = QueryCache::new();
    let is_valid = v.validate_w_cache(&mut cache);
    assert!(is_valid)
}

trait ValidateAgainst<OtherT> {
    fn validate_against(&self, other: &OtherT) -> bool;
}

impl<OtherT, T> ValidateAgainst<OtherT> for T
where
    T: Pufferfish,
    for<'a> OtherT: PuffQuery<'a>,
    // for<'b> <T as PuffQuery<'a>>::HitsT: DecodeHit<'b>,
    // for<'b> <OtherT as PuffQuery<'o>>::HitsT: DecodeHit<'b>,
{
    // impl<T: DebugGetRefPos + PufferfishBase> ValidateAgainst<T> for SparseContigDenseIndex {
    fn validate_against(&self, other: &OtherT) -> bool {
        // Only check if self has reference sequence encoded
        if !self.has_ref_seq() {
            // TODO message to tell user to get set pufferfish flag for ref seq
            panic!("Cannot validate index without stored reference sequence");
        }

        for (_ref_id, refseq) in self.iter_refs().enumerate() {
            for (_pos, kmer) in refseq.iter_kmers().enumerate() {
                let canon_kmer = CanonicalKmer::from(kmer);
                let hits = self.get_ref_pos(&canon_kmer).unwrap().as_vec();
                let other_hits = other.get_ref_pos(&canon_kmer).unwrap().as_vec();

                if hits != other_hits {
                    panic!("WARN: validation failed");
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod test_cached {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    // use crate::index::sparse_contig_index::SparseDenseIndex;
    use crate::index::{DenseIndex, SparseIndex};
    use crate::test_utils::*;

    // TODO validate cached for WMSparseSparse and WMSparseDense
    #[test]
    fn small_txome_dense() {
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        // pi.validate();
        assert_is_valid_cached(&pi);
    }

    #[test]
    fn small_txome_sparse() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        // pi.validate();
        assert_is_valid_cached(&pi);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::index::{DenseIndex, SparseIndex};
    use crate::test_utils::*;

    #[test]
    fn validate_tiny_index() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        // assert!(pi.validate());
        assert_is_valid(&pi);
    }

    #[test]
    fn validate_tiny_refs_index() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        // assert!(pi.validate());
        assert_is_valid(&pi);
    }

    #[test]
    fn validate_yeast_index() {
        let p = to_abs_path(YEAST_CHR01_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        // assert!(pi.validate());
        assert_is_valid(&pi);
    }

    #[test]
    fn small_txome_dense() {
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        // pi.validate();
        assert_is_valid(&pi);
    }

    #[test]
    fn small_txome_sparse() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        // pi.validate();
        assert_is_valid(&pi);
    }
}

#[cfg(test)]
mod validate_against_tests {
    use super::*;
    use crate::index::{DenseIndex, SparseIndex};
    // use crate::index::sparse_contig_index::SparseContigDenseIndex;
    use crate::cpp::DeserializeFromCpp;
    use crate::test_utils::*;
    #[test]
    fn sparse_vs_dense() {
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let dpi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let spi = SparseIndex::deserialize_from_cpp(&p).unwrap();

        assert!(dpi.validate_against(&spi));
        assert!(spi.validate_against(&dpi));
        assert!(dpi.validate_against(&dpi));
        assert!(spi.validate_against(&spi));
    }
}
