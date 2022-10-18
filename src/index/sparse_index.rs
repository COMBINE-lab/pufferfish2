use super::base_index::BaseIndex;
use super::compact::FromCompact;

use super::projected_hits::ProjectedHits;
use super::refseq_iter::RefIterator;
use super::{
    CachedQuery, DenseContigTable, Info, MappedOrientation, PuffQuery, PufferfishBase,
    PufferfishFilePaths, PufferfishType, QueryCache,
};
use crate::boophf::MPHF;
use crate::cpp::DeserializeFromCpp;
use kmers::naive_impl::{CanonicalKmer, MatchType};

use delegate::delegate;
use simple_sds::bit_vector::BitVector;
use simple_sds::int_vector::IntVector;
use simple_sds::raw_vector::AccessRaw;

use simple_sds::ops::{Access, BitVec, Rank, Select};
use std::io::Result;
use std::path::Path;

#[derive(Clone)]
pub struct SparseIndex {
    // all fields are pub(super) for easy unpacking in super/crate::index
    pub(super) base: BaseIndex,
    pub(super) ctg_table: DenseContigTable,
    pub(super) sampled_vec: BitVector,   // is position sampled
    pub(super) canonical_vec: BitVector, // orientation for non-sampled kmers
    pub(super) direction_vec: BitVector, // which direction to walk
    pub(super) ext_sizes: IntVector,     // extension sizes
    pub(super) ext_bases: IntVector,     // the extension nucleotides
    pub(super) sampled_pos: IntVector,
    pub(super) sample_size: usize,
    pub(super) extension_size: usize,
}

impl DeserializeFromCpp for SparseIndex {
    fn deserialize_from_cpp<P: AsRef<Path>>(dpath: P) -> Result<Self> {
        if let Ok(base_struct) = BaseIndex::deserialize_from_cpp(&dpath) {
            let files = PufferfishFilePaths::new(&dpath);

            let info = Info::load(files.info_json);
            assert_eq!(info.sampling_type, PufferfishType::Sparse); //TODO return an error instaed

            let sampled_pos = IntVector::from_compact_serialized(files.sample_pos)?;
            let canonical_vec = BitVector::from_compact_serialized(files.canonical)?;
            let direction_vec = BitVector::from_compact_serialized(files.direction)?;
            let ext_sizes = IntVector::from_compact_serialized(files.extension_lengths)?;
            let ext_bases = IntVector::from_compact_serialized(files.extension_bases)?;
            let sampled_vec = {
                let mut sv = BitVector::from_compact_serialized(files.presence)?;
                sv.enable_rank();
                sv.enable_select();
                sv
            };

            let sample_size = info.sample_size.unwrap();
            let extension_size = info.extension_size.unwrap();

            let ctg_table = DenseContigTable::deserialize_from_cpp(&dpath)?;

            Ok(Self {
                base: base_struct,
                ctg_table,
                sampled_vec,
                canonical_vec,
                direction_vec,
                ext_sizes,
                ext_bases,
                sampled_pos,
                sample_size,
                extension_size,
            })
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "could not load index.",
            ))
        }
    }
}
impl SparseIndex {
    // pub fn get_contig_kmer_u64(&self, contig_id: usize, offset: usize) -> u64 {
    //     self.base.get_contig_kmer_u64(contig_id, offset)
    // }

    pub fn get_ctg_table(&self) -> &DenseContigTable {
        &self.ctg_table
    }

    #[allow(dead_code)]
    #[inline]
    fn get_sample_size(&self) -> usize {
        self.sample_size
    }

    fn get_ref_pos_helper(
        &self,
        mer: &CanonicalKmer,
        pos: usize,
        did_walk: bool,
    ) -> Option<ProjectedHits> {
        if pos <= self.base.last_seq_pos {
            let fk = self.base.seq.get_kmer_u64(pos, self.base.k as u8);
            let keq = mer.get_word_equivalency(fk);

            if keq == MatchType::NoMatch {
                return None;
            }

            let rank_interval = if did_walk {
                unsafe { self.base.bv.as_ref().int(pos, (self.base.k as usize) - 1) }
            } else {
                0
            };

            if rank_interval > 0 {
                return None;
            }

            let contig_id = self.base.bv.rank(pos) as usize;

            let offset = self.contig_start_pos(contig_id);
            let contig_len = self.contig_len(contig_id);

            let mapped_orientation = match keq {
                MatchType::IdentityMatch => MappedOrientation::Forward,
                MatchType::TwinMatch => MappedOrientation::Backward,
                // _ => MappedOrientation::Forward,
                _ => panic!("Cannot convert MatchType::NoMatch to MappedOrientation"),
            };

            let offset = pos - offset;
            let hits = self.ctg_table.get_encoded_contig_occs(contig_id);
            Some(ProjectedHits {
                kmer: mer.get_fw_mer(),
                mapped_orientation,
                offset,
                hits,
                contig_len,
                contig_id,
            })
        } else {
            // trace!("exiting helper {} {}", pos, self.base.last_seq_pos);
            None
        }
    }

    fn get_ref_pos_helper_with_cache(
        &self,
        mer: &CanonicalKmer,
        pos: usize,
        did_walk: bool,
        qc: &mut QueryCache,
    ) -> Option<ProjectedHits> {
        if pos <= self.base.last_seq_pos {
            let fk = self.base.seq.get_kmer_u64(pos, self.base.k as u8);
            let keq = mer.get_word_equivalency(fk);

            if keq == MatchType::NoMatch {
                // trace!("no match by equivalency");
                return None;
            }

            let rank_interval = if did_walk {
                unsafe { self.base.bv.as_ref().int(pos, (self.base.k as usize) - 1) }
            } else {
                0
            };

            if rank_interval > 0 {
                // trace!("no match by rank_interval");
                return None;
            }

            let mapped_orientation = match keq {
                MatchType::IdentityMatch => MappedOrientation::Forward,
                MatchType::TwinMatch => MappedOrientation::Backward,
                _ => MappedOrientation::Forward,
            };

            let contig_id;
            let offset;
            let contig_len;

            let same_contig_as_cache =
                (qc.contig_start <= pos) && (pos < qc.contig_start + qc.contig_len);

            if same_contig_as_cache {
                contig_id = qc.contig_id;
                contig_len = qc.contig_len;
                offset = pos - qc.contig_start;
            } else {
                contig_id = self.base.bv.rank(pos) as usize;

                // start position of contig on seq
                let ctg_start_pos = if contig_id == 0 {
                    0
                } else {
                    self.base.bv.select(contig_id - 1).unwrap() as usize + 1
                };

                contig_len = self.base.bv.select(contig_id).unwrap() + 1 - ctg_start_pos;
                offset = pos - ctg_start_pos;

                qc.contig_len = contig_len;
                qc.contig_id = contig_id;
                qc.contig_start = ctg_start_pos;
                // dbg!("cache miss");
            }

            let hits = self.ctg_table.get_encoded_contig_occs(contig_id);
            Some(ProjectedHits {
                kmer: mer.get_fw_mer(),
                mapped_orientation,
                offset,
                hits,
                contig_len,
                contig_id,
            })
        } else {
            None
        }
    }
}

impl PufferfishBase for SparseIndex {
    fn num_total_occs(&self) -> usize {
        self.ctg_table.num_total_occs()
    }
    fn num_ctg_occs(&self, contig_id: usize) -> usize {
        self.ctg_table.num_ctg_occs(contig_id)
    }

    fn get_base_index(&self) -> &BaseIndex {
        &self.base
    }

    delegate! {
        to self.ctg_table {
            fn get_ref_name(&self, id: usize) -> Option<&String>;
            fn get_ref_names(&self) -> &[String];
            fn get_encoded_contig_occs(&self, contig_id: usize) -> &[u64];
        }
    }
    fn to_info(&self) -> Info {
        let base = &self.base;
        Info {
            sampling_type: PufferfishType::SparseRS,
            index_version: base.index_version,
            reference_gfa: base.reference_gfa.clone(),
            kmer_size: base.k,
            num_kmers: base.num_kmers(),
            num_contigs: base.num_contigs(),
            seq_len: base.seq_len(),
            have_ref_seq: base.ref_seq.is_some(),
            have_edge_vec: base.have_edge_vec,
            seq_hash: base.seq_hash.clone(),
            name_hash: base.name_hash.clone(),
            seq_hash_512: base.seq_hash_512.clone(),
            name_hash_512: base.name_hash_512.clone(),
            decoy_seq_hash: base.decoy_seq_hash.clone(),
            decoy_name_hash: base.decoy_name_hash.clone(),
            num_decoys: base.num_decoys,
            sample_size: Some(self.sample_size), // Putting into enum would break c++ compatiblity)
            extension_size: Some(self.extension_size),
            first_decoy_index: base.first_decoy_index,
            keep_duplicates: base.keep_duplicates,
        }
    }

    fn iter_refs(&self) -> RefIterator<Self> {
        // fn iter_refs<'a: 'b, 'b>(&'a self) -> RefIterator<'b, Self> {
        if !self.has_ref_seq() {
            // TODO: We may want to move this panic inside of the kmer iterators
            // if the BaseIndexReference supports more functionalities that don't
            // require the sequence
            panic!("Cannot iterate over references with no stored refseq")
        }
        RefIterator::new(self)
    }
}

impl<'a> PuffQuery<'a> for SparseIndex {
    type HitsT = ProjectedHits<'a>;
    fn get_ref_pos(&'a self, kmer: &CanonicalKmer) -> Option<Self::HitsT> {
        let mut did_walk = false;

        let kw = kmer.get_canonical_word();
        let mut km = kmer.clone();
        if !km.is_fw_canonical() {
            km.swap();
        }

        let idx = self.base.mphf.lookup(&kw)? as usize;

        let pos: usize;
        if self.sampled_vec.get(idx) {
            pos = self.sampled_pos.get(self.sampled_vec.rank(idx)) as usize;
            // trace!("is sampled w pos {}", pos);
        } else {
            // trace!("is not sampled");
            did_walk = true;
            let mut signed_shift = 0i64;
            //let mut in_loop = 0;

            let current_rank = self.sampled_vec.rank(idx);
            let extension_pos = idx - current_rank;
            let extension_word = self.ext_bases.get(extension_pos);

            if !self.canonical_vec.get(extension_pos) && km.is_fw_canonical() {
                km.swap();
            }

            let shift_fw = self.direction_vec.get(extension_pos);
            let llimit =
                self.extension_size - (self.ext_sizes.get(extension_pos) as usize + 1usize);

            if shift_fw {
                for i in (llimit + 1..=self.extension_size).rev() {
                    let ssize = 2 * (i - 1);
                    let curr_code = (extension_word & (0x3 << ssize)) >> ssize;
                    km.append_base(curr_code);
                    signed_shift -= 1;
                }
            } else {
                for i in (llimit + 1..=self.extension_size).rev() {
                    let ssize = 2 * (i - 1);
                    let curr_code = (extension_word & (0x3 << ssize)) >> ssize;
                    km.prepend_base(curr_code);
                    signed_shift += 1;
                }
            }

            let kw = km.get_canonical_word();
            let idx = self.base.mphf.lookup(&kw)? as usize;
            if !self.sampled_vec.get(idx) {
                return None;
            }

            let current_rank = self.sampled_vec.rank(idx);
            let sample_pos = self.sampled_pos.get(current_rank);

            pos = ((sample_pos as i64) + signed_shift) as usize;
        }

        self.get_ref_pos_helper(kmer, pos, did_walk)
    }
}

impl<'a> CachedQuery<'a, QueryCache> for SparseIndex {
    type HitsT = ProjectedHits<'a>;
    fn get_ref_pos_with_cache(
        &'a self,
        kmer: &CanonicalKmer,
        qc: &mut QueryCache,
    ) -> Option<Self::HitsT> {
        let mut did_walk = false;

        let kw = kmer.get_canonical_word();
        let mut km = kmer.clone();
        if !km.is_fw_canonical() {
            km.swap();
        }

        // idiomatic to use ? here; see https://rust-lang.github.io/rust-clippy/master/index.html#question_mark
        let idx = self.base.mphf.lookup(&kw)? as usize;

        let pos: usize;
        if self.sampled_vec.get(idx) {
            pos = self.sampled_pos.get(self.sampled_vec.rank(idx)) as usize;
        } else {
            did_walk = true;
            let mut signed_shift = 0i64;
            //let mut in_loop = 0;

            let current_rank = self.sampled_vec.rank(idx);
            let extension_pos = idx - current_rank;
            let extension_word = self.ext_bases.get(extension_pos);

            if !self.canonical_vec.get(extension_pos) && km.is_fw_canonical() {
                km.swap();
            }

            let shift_fw = self.direction_vec.get(extension_pos);
            let llimit =
                self.extension_size - (self.ext_sizes.get(extension_pos) as usize + 1usize);

            if shift_fw {
                for i in (llimit + 1..=self.extension_size).rev() {
                    let ssize = 2 * (i - 1);
                    let curr_code = (extension_word & (0x3 << ssize)) >> ssize;
                    km.append_base(curr_code);
                    signed_shift -= 1;
                }
            } else {
                for i in (llimit + 1..=self.extension_size).rev() {
                    let ssize = 2 * (i - 1);
                    let curr_code = (extension_word & (0x3 << ssize)) >> ssize;
                    km.prepend_base(curr_code);
                    signed_shift += 1;
                }
            }

            let kw = km.get_canonical_word();
            let idx = self.base.mphf.lookup(&kw)? as usize;
            if !self.sampled_vec.get(idx) {
                return None;
            }

            let current_rank = self.sampled_vec.rank(idx);
            //in_loop += 1;
            let sample_pos = self.sampled_pos.get(current_rank);

            pos = ((sample_pos as i64) + signed_shift) as usize;
        }

        self.get_ref_pos_helper_with_cache(kmer, pos, did_walk, qc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::PufferfishType;
    use crate::test_utils::*;

    #[test]
    fn pufferfish_type() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        assert_eq!(pi.to_info().sampling_type, PufferfishType::SparseRS);
    }
}
