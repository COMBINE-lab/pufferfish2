use bincode;

use log::debug;
use simple_sds::bit_vector::BitVector;
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Access, BitVec, Rank};
use simple_sds::raw_vector::AccessRaw;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

use super::contig_samplers::ContigSamplingStrategy;
use super::occs::EncodedOccs;
use super::wm_sparse_contig_table::{SparseContigTable, SparseContigTableBuilder};
use super::wm_sparse_projected_hits::SparseProjHits;
use super::StreamingCache;
use super::DEFAULT_CONTIG_SAMPLIG;

use crate::boophf::MPHF;
use crate::index::base_index::BaseIndex;
use crate::index::info::Info;
use crate::index::refseq_iter::RefIterator;
use crate::index::{MappedOrientation, PuffQuery, PufferfishBase, PufferfishType, SparseIndex};
use crate::io::{DeserializeFrom, SerializeTo};
use crate::serde_ext::{AsSerialize, FromDeserialize};
use kmers::naive_impl::{CanonicalKmer, MatchType};

impl<'a> SparseSparseIndex {
    pub fn grp_w_cache(
        &'a self,
        kmer: &CanonicalKmer,
        cache: &mut StreamingCache,
    ) -> Option<SparseProjHits<'a, Self>> {
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
            let sample_pos = self.sampled_pos.get(current_rank);
            pos = ((sample_pos as i64) + signed_shift) as usize;
        }

        self.grp_w_cache_helper(kmer, pos, did_walk, cache)
    }

    fn grp_w_cache_helper(
        &'a self,
        mer: &CanonicalKmer,
        pos: usize,
        did_walk: bool,
        cache: &mut StreamingCache,
    ) -> Option<SparseProjHits<'a, Self>> {
        if pos > self.base.last_seq_pos {
            return None;
        }

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

        let mapped_orientation = match keq {
            MatchType::IdentityMatch => MappedOrientation::Forward,
            MatchType::TwinMatch => MappedOrientation::Backward,
            // _ => MappedOrientation::Forward,
            _ => panic!("Cannot convert MatchType::NoMatch to MappedOrientation"),
        };

        let contig_id;
        let contig_len;
        let offset;
        let occs;

        let same_contig_as_cache =
            (cache.contig_start <= pos) && (pos < cache.contig_start + cache.contig_len);

        if same_contig_as_cache {
            // Cache hit
            // debug!("HIT CONTIG ID");
            contig_id = cache.contig_id;
            contig_len = cache.contig_len;
            offset = pos - cache.contig_start;
            if let Some(inner_occs) = &cache.contig_occs {
                // debug!("HIT DECODED_OCCS");
                // Cache hit *and* cache knows encoded occs as vec of mapped_ref_pos.
                occs = EncodedOccs::NotEncoded(inner_occs.clone());
            } else {
                occs = self.ctg_table.get_encoded_contig_occs(contig_id);
            }
        } else {
            // Cache miss.
            contig_id = self.base.bv.rank(pos) as usize;

            let start_pos = self.contig_start_pos(contig_id);
            contig_len = self.contig_len(contig_id);
            offset = pos - start_pos;
            occs = self.ctg_table.get_encoded_contig_occs(contig_id);

            cache.contig_id = contig_id;
            cache.contig_start = start_pos;
            cache.contig_len = contig_len;
            cache.contig_occs = None;
            // None
        }
        let hits = SparseProjHits {
            index: self,
            kmer: mer.get_fw_mer(),
            mapped_orientation,
            offset,
            hits: occs,
            contig_len,
            contig_id,
        };
        Some(hits)
        // } else {
        //     None
        // }
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct SparseSparseIndex {
    base: BaseIndex,
    ctg_table: SparseContigTable,

    // TODO pub crate?
    pub contig_sampling_strategy: ContigSamplingStrategy,
    sampled_vec: BitVector,   // is position sampled
    canonical_vec: BitVector, // orientation for non-sampled kmers
    direction_vec: BitVector, // which direction to walk
    ext_sizes: IntVector,     // extension sizes
    ext_bases: IntVector,     // the extension nucleotides
    sampled_pos: IntVector,

    pub sampled_wm_thresh: usize, // Threshhold # of occs to get WM support
    pub nonsamp_wm_thresh: usize, // Threshhold # of occs to get WM support

    pub sample_size: usize,
    pub extension_size: usize,
}

impl PufferfishBase for SparseSparseIndex {
    fn num_total_occs(&self) -> usize {
        self.ctg_table.num_total_occs()
    }
    fn num_ctg_occs(&self, contig_id: usize) -> usize {
        self.ctg_table.num_ctg_occs(contig_id)
    }

    fn get_encoded_contig_occs(&self, _contig_id: usize) -> &[u64] {
        panic!("Not available for WMSparseSparseIndex")
    }

    fn get_base_index(&self) -> &BaseIndex {
        &self.base
    }

    fn to_info(&self) -> Info {
        let base = &self.base;
        Info {
            sampling_type: PufferfishType::SparseSparseV2(self.contig_sampling_strategy),
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

    fn get_ref_name(&self, _ref_id: usize) -> Option<&String> {
        todo!()
    }
    fn get_ref_names(&self) -> &[String] {
        todo!()
    }
    // fn get_ref_name(&self) -> String { }

    fn iter_refs(&self) -> RefIterator<Self> {
        if !self.has_ref_seq() {
            // TODO: We may want to move this panic inside of the kmer iterators
            // if the BaseIndexReference supports more functionalities that don't
            // require the sequence
            panic!("Cannot iterate over references with no stored refseq")
        }
        RefIterator::new(self)
    }
}

impl<'a> PuffQuery<'a> for SparseSparseIndex {
    type HitsT = SparseProjHits<'a, Self>;

    fn get_ref_pos(&self, kmer: &CanonicalKmer) -> Option<SparseProjHits<Self>> {
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
            let sample_pos = self.sampled_pos.get(current_rank);
            pos = ((sample_pos as i64) + signed_shift) as usize;
        }

        self.get_ref_pos_helper(kmer, pos, did_walk)
    }
}

impl SparseSparseIndex {
    pub fn contig_is_sampled(&self, contig_id: usize) -> bool {
        self.ctg_table.is_sampled(contig_id)
    }

    pub fn contig_num_occs(&self, contig_id: usize) -> usize {
        self.ctg_table.num_ctg_occs(contig_id)
    }

    fn get_ref_pos_helper(
        &self,
        mer: &CanonicalKmer,
        pos: usize,
        did_walk: bool,
    ) -> Option<SparseProjHits<Self>> {
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

            // // start position of contig on seq
            // let offset = if contig_id == 0 {
            //     0
            // } else {
            //     self.base.bv.select(contig_id - 1).unwrap() as usize + 1
            // };

            // let contig_len = self.base.bv.select(contig_id).unwrap() + 1 - offset;

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
            // let kth_bases = if let EncodedOccs::NotSampledWM(_) = hits {
            //     None
            // } else {
            //     None
            // };
            Some(SparseProjHits {
                index: self,
                kmer: mer.get_fw_mer(),
                mapped_orientation,
                offset,
                hits,
                contig_len,
                contig_id,
                // kth_bases,
            })
        } else {
            None
        }
    }

    pub fn num_sampled_contigs(&self) -> usize {
        self.ctg_table.num_sampled()
    }

    pub fn num_not_sampled_contigs(&self) -> usize {
        self.ctg_table.num_not_sampled()
    }
}

impl SparseSparseIndex {
    pub fn from_sparse_index(
        pi: SparseIndex,
        pop_sampled_thresh: usize,
        pop_nonsamp_thresh: usize,
    ) -> Self {
        Self::from_sparse_index_w_strategy(
            pi,
            DEFAULT_CONTIG_SAMPLIG,
            pop_sampled_thresh,
            pop_nonsamp_thresh,
        )
    }

    pub fn from_sparse_index_w_strategy(
        pi: SparseIndex,
        strategy: ContigSamplingStrategy,
        sampled_wm_thresh: usize,
        nonsamp_wm_thresh: usize,
    ) -> Self {
        assert!(pi.base.ref_seq.is_some());

        let ctg_table =
            SparseContigTableBuilder::new(&pi, strategy, sampled_wm_thresh, nonsamp_wm_thresh)
                .build();
        Self {
            ctg_table,
            base: pi.base,
            contig_sampling_strategy: strategy,

            sampled_wm_thresh,
            nonsamp_wm_thresh,

            sampled_vec: pi.sampled_vec,
            canonical_vec: pi.canonical_vec,
            direction_vec: pi.direction_vec,
            ext_sizes: pi.ext_sizes,
            ext_bases: pi.ext_bases,
            sampled_pos: pi.sampled_pos,
            sample_size: pi.sample_size,
            extension_size: pi.extension_size,
        }
    }

    // Serialization

    pub fn deserialize_from_parts<A: AsRef<Path>, B: AsRef<Path>, C: AsRef<Path>>(
        base_dir: A,
        k2u_dir: B,
        ctable_dir: C,
    ) -> bincode::Result<Self> {
        // Expects
        // Serializes to dirs like
        // dir
        //  |- mphf.bin
        //  |- base_index.bin
        //  |- info.json
        //  |- {ctable_dir}/
        //       |- ... ctable files
        //  |- {k2u_dir}/...
        //       |- ... k2u files
        let dir = base_dir.as_ref();
        let k2u_dir = k2u_dir.as_ref();
        let ctable_dir = ctable_dir.as_ref();
        let fps = SparseSparseIndexPaths::new_from_parts(dir, k2u_dir, ctable_dir);
        Self::deserialize_from_paths(&fps)
    }

    fn deserialize_from_paths(fps: &SparseSparseIndexPaths) -> bincode::Result<Self> {
        debug!("Loading info");
        // let base_info = BaseInfo::load(&fps.info);
        let ctable_info = CTableInfo::load(&fps.ctable_info);
        let k2u_info = K2UInfo::load(&fps.k2u_info);

        // Load and check contig sampling strategy
        let contig_sampling_strategy;
        let nonsamp_wm_thresh;
        let sampled_wm_thresh;
        if let CTableInfo::SparseWM {
            sampled_wm_thresh: sampled_t,
            nonsamp_wm_thresh: nonsamp_t,
            sampling,
        } = ctable_info
        {
            contig_sampling_strategy = sampling;
            nonsamp_wm_thresh = nonsamp_t;
            sampled_wm_thresh = sampled_t;
        } else {
            panic!(
                "Cannot deserialize WMSparseSparseIndex with {:?}",
                ctable_info
            )
        };

        // Load and check pos sampling strategy
        let extension_size;
        let sample_size;
        if let K2UInfo::Sparse {
            sample_size: s,
            extension_size: e,
        } = k2u_info
        {
            sample_size = s;
            extension_size = e;
        } else {
            panic!("Cannot deserialize WMSparseSparseIndex with {:?}", k2u_info)
        };

        debug!("Loading sparse contig table");
        let ctg_table = SparseContigTable::deserialize_from(&fps.ctable_dir)?;

        // let extension_size = info.extension_size.unwrap();
        // let _sample_size = info.sample_size.unwrap();

        debug!("Loading base index");
        let f = File::open(&fps.base)?;
        let r = BufReader::new(f);
        let mut base: BaseIndex = bincode::deserialize_from(r)?;

        debug!("Loading MPHF");
        let f = File::open(&fps.mphf)?;
        let r = BufReader::new(f);
        let mphf = bincode::deserialize_from(r)?;
        base.mphf = mphf;

        debug!("Loading sample pos");
        let f = File::open(&fps.sampled_pos)?;
        let r = BufReader::new(f);
        let de = bincode::deserialize_from(r)?;
        let sampled_pos = IntVector::from_deserialized(de);

        debug!("Loading sample vec");
        let f = File::open(&fps.sampled_vec)?;
        let r = BufReader::new(f);
        let de = bincode::deserialize_from(r)?;
        let sampled_vec = BitVector::from_deserialized(de);

        debug!("Loading ext sizes");
        let f = File::open(&fps.ext_sizes)?;
        let r = BufReader::new(f);
        let de = bincode::deserialize_from(r)?;
        let ext_sizes = IntVector::from_deserialized(de);

        debug!("Loading ext bases");
        let f = File::open(&fps.ext_bases)?;
        let r = BufReader::new(f);
        let de = bincode::deserialize_from(r)?;
        let ext_bases = IntVector::from_deserialized(de);

        debug!("Loading canonical vec");
        let f = File::open(&fps.canonical_vec)?;
        let r = BufReader::new(f);
        let de = bincode::deserialize_from(r)?;
        let canonical_vec = BitVector::from_deserialized(de);

        debug!("Loading direction vec");
        let f = File::open(&fps.direction_vec)?;
        let r = BufReader::new(f);
        let de = bincode::deserialize_from(r)?;
        let direction_vec = BitVector::from_deserialized(de);

        let index = Self {
            base,
            ctg_table,

            contig_sampling_strategy,
            sampled_wm_thresh,
            nonsamp_wm_thresh,

            sampled_pos,
            sampled_vec,
            ext_bases,
            ext_sizes,
            canonical_vec,
            direction_vec,

            sample_size,
            extension_size,
        };
        Ok(index)
    }

    pub fn lazy_serialize_to<P: AsRef<Path>>(&self, dir: P) -> bincode::Result<()> {
        // Serializes to dirs (with auto ctable_dir and k2u_dir names)
        // dir
        //  |- {mphf, base_index}.bin
        //  |- info.json
        //  |- {ctable_dir}
        //       |- ... ctable files
        //  |- {k2u_dir}/...
        //       |- ... k2u files

        let v2_info = V2Info::from_wm_sparse_sparse_index(self);

        let dir = dir.as_ref();
        let fps = SparseSparseIndexPaths::new_from_parts(
            dir,
            dir.join(v2_info.k2u_info.to_dir_name()),
            dir.join(v2_info.ctable_info.to_dir_name()),
        );

        self.serialize_base_lazy(fps.clone())?;
        self.serialize_ctable_lazy(fps.clone())?;
        self.serialize_k2u_lazy(fps)?;
        Ok(())
    }

    fn serialize_base_lazy(&self, fps: SparseSparseIndexPaths) -> bincode::Result<()> {
        let v2_info = V2Info::from_wm_sparse_sparse_index(self).base_info;

        if fps.info.is_file() {
            let old_info = BaseInfo::load(&fps.info);
            debug!(
                "Skipping serialization of BaseIndex, MPHF, and Info, {:?} found on disk",
                fps.info
            );

            assert_eq!(
                old_info, v2_info,
                "Cannot skip serialization, BaseInfo mismatch -- something went wrong"
            );

            assert!(fps.base.is_file(), "Base missing");
            assert!(fps.mphf.is_file(), "MPHF missing");
        } else {
            debug!("Serializing BaseIndex, MPHF, and Info");
            assert!(!fps.base.is_file(), "Base already exists");
            assert!(!fps.mphf.is_file(), "MPHF already exists");
            std::fs::create_dir_all(fps.dir)?;
            // INFO
            let f = File::create(fps.info)?;
            let w = BufWriter::new(f);
            // serde_json::to_writer_pretty(w, &self.to_info()).unwrap();
            serde_json::to_writer_pretty(w, &v2_info).unwrap();

            // BASE + MPHF
            let f = File::create(fps.base)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.base)?;

            let f = File::create(fps.mphf)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.base.mphf)?;
        }

        Ok(())
    }
    fn serialize_k2u_lazy(&self, fps: SparseSparseIndexPaths) -> bincode::Result<()> {
        let v2_info = V2Info::from_wm_sparse_sparse_index(self);
        if fps.k2u_info.is_file() {
            let old_info = K2UInfo::load(&fps.k2u_info);
            debug!(
                "Skipping serialization of K2U data structures, {:?} found on disk",
                fps.k2u_info
            );
            assert_eq!(
                old_info, v2_info.k2u_info,
                "Cannot skip serialization, K2U type mismatch -- something went wrong"
            );
            assert!(fps.sampled_pos.is_file());
            assert!(fps.sampled_vec.is_file());
            assert!(fps.canonical_vec.is_file());
            assert!(fps.direction_vec.is_file());
            assert!(fps.ext_sizes.is_file());
            assert!(fps.ext_bases.is_file());
        } else {
            debug!("Serializing K2U data structures");
            assert!(!fps.sampled_pos.is_file());
            assert!(!fps.sampled_vec.is_file());
            assert!(!fps.canonical_vec.is_file());
            assert!(!fps.direction_vec.is_file());
            assert!(!fps.ext_sizes.is_file());
            assert!(!fps.ext_bases.is_file());
            std::fs::create_dir_all(fps.k2u_dir)?;

            let f: File = File::create(fps.k2u_info)?;
            serde_json::to_writer_pretty(f, &v2_info.k2u_info).unwrap();

            let f = File::create(fps.sampled_vec)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.sampled_vec.as_serialize())?;

            let f = File::create(fps.sampled_pos)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.sampled_pos.as_serialize())?;

            let f = File::create(fps.canonical_vec)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.canonical_vec.as_serialize())?;

            let f = File::create(fps.direction_vec)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.direction_vec.as_serialize())?;

            let f = File::create(fps.ext_sizes)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.ext_sizes.as_serialize())?;

            let f = File::create(fps.ext_bases)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.ext_bases.as_serialize())?;
        }
        Ok(())
    }

    fn serialize_ctable_lazy(&self, fps: SparseSparseIndexPaths) -> bincode::Result<()> {
        let v2_info = V2Info::from_wm_sparse_sparse_index(self);
        if fps.ctable_info.is_file() {
            let old_info = CTableInfo::load(&fps.ctable_info);
            debug!(
                "Skipping serialization of contig table data structures, {:?} found on disk",
                fps.ctable_info
            );
            assert_eq!(
                old_info, v2_info.ctable_info,
                "Cannot skip serialization, contig table type mismatch -- something went wrong"
            );
            log::warn!("No checks for sparse ctable files are implemented yet!");
            assert!(fps.ctable_dir.is_dir());
        } else {
            debug!("Serializing contig table");
            assert!(!fps.ctable_dir.exists());

            std::fs::create_dir_all(&fps.ctable_dir)?;

            let f: File = File::create(fps.ctable_info)?;
            serde_json::to_writer_pretty(f, &v2_info.ctable_info).unwrap();

            self.ctg_table.serialize_to(&fps.ctable_dir)?;
        }
        Ok(())
    }
}

use crate::index::info::{BaseInfo, CTableInfo, K2UInfo, V2Info};

// Serialization
impl<P: AsRef<Path>> SerializeTo<P> for SparseSparseIndex {
    fn serialize_to(&self, dir: P) -> bincode::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir(dir)?;
        let fps = SparseSparseIndexPaths::new(dir);

        let v2_info = V2Info::from_wm_sparse_sparse_index(self);
        let f = File::create(fps.info)?;
        serde_json::to_writer_pretty(f, &v2_info.base_info).unwrap();
        let f = File::create(fps.k2u_info)?;
        serde_json::to_writer_pretty(f, &v2_info.k2u_info).unwrap();
        let f = File::create(fps.ctable_info)?;
        serde_json::to_writer_pretty(f, &v2_info.ctable_info).unwrap();

        self.ctg_table.serialize_to(&dir)?;

        let f = File::create(fps.base)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.base)?;

        let f = File::create(fps.mphf)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.base.mphf)?;

        let f = File::create(fps.sampled_vec)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.sampled_vec.as_serialize())?;

        let f = File::create(fps.sampled_pos)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.sampled_pos.as_serialize())?;

        let f = File::create(fps.canonical_vec)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.canonical_vec.as_serialize())?;

        let f = File::create(fps.direction_vec)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.direction_vec.as_serialize())?;

        let f = File::create(fps.ext_sizes)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.ext_sizes.as_serialize())?;

        let f = File::create(fps.ext_bases)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.ext_bases.as_serialize())?;
        Ok(())
    }
}

impl<P: AsRef<Path>> DeserializeFrom<P> for SparseSparseIndex {
    fn deserialize_from(dir: P) -> bincode::Result<Self> {
        let dir = dir.as_ref();
        let fps = SparseSparseIndexPaths::new(dir);
        Self::deserialize_from_paths(&fps)
    }
}

#[derive(Clone)]
struct SparseSparseIndexPaths {
    dir: PathBuf,
    k2u_dir: PathBuf,
    ctable_dir: PathBuf,

    info: PathBuf,
    mphf: PathBuf,
    base: PathBuf,

    k2u_info: PathBuf,
    ctable_info: PathBuf,

    sampled_vec: PathBuf,   // is position sampled
    canonical_vec: PathBuf, // orientation for non-sampled kmers
    direction_vec: PathBuf, // which direction to walk
    ext_sizes: PathBuf,
    ext_bases: PathBuf, // the extension nucleotides
    sampled_pos: PathBuf,
}

impl SparseSparseIndexPaths {
    fn new_from_parts<P1: AsRef<Path>, P2: AsRef<Path>, P3: AsRef<Path>>(
        dir: P1,
        k2u_dir: P2,
        ctable_dir: P3,
    ) -> Self {
        let dir = dir.as_ref();
        let k2u_dir = k2u_dir.as_ref();
        let ctable_dir = ctable_dir.as_ref();

        Self {
            dir: dir.to_path_buf(),
            ctable_dir: ctable_dir.to_path_buf(),
            k2u_dir: k2u_dir.to_path_buf(),

            info: dir.join("info.json"),
            k2u_info: k2u_dir.join("k2u_info.json"),
            ctable_info: ctable_dir.join("ctable_info.json"),

            mphf: dir.join("mphf.bin"),
            base: dir.join("base_index.bin"),

            sampled_pos: k2u_dir.join("sample_pos.bin"),
            canonical_vec: k2u_dir.join("canonical.bin"),
            direction_vec: k2u_dir.join("direction.bin"),
            ext_sizes: k2u_dir.join("extensionSize.bin"),
            ext_bases: k2u_dir.join("extension.bin"),
            sampled_vec: k2u_dir.join("presence.bin"),
        }
    }

    fn new<P: AsRef<Path>>(dir: P) -> Self {
        let dir = dir.as_ref();
        Self::new_from_parts(dir, dir, dir)
    }
}

#[cfg(test)]
mod serialization_tests {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::io::{DeserializeFrom, SerializeTo};
    use crate::test_utils::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn ser_de_small() {
        let dir = temp_file_name("");
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseSparseIndex::from_sparse_index(pi, 0, 0);

        spi.serialize_to(&dir).unwrap();

        let spi_loaded = SparseSparseIndex::deserialize_from(&dir).unwrap();

        std::fs::remove_dir_all(dir).unwrap();

        assert_eq!(spi, spi_loaded);
    }
}

#[cfg(test)]
#[cfg(test)]
mod wm_tests {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::index::sparse_contig_index::DEFAULT_CONTIG_SAMPLIG;
    use crate::index::validate::*;
    use crate::test_utils::*;
    #[test]
    fn pufferfish_type() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseSparseIndex::from_sparse_index(pi, 0, 0);
        assert_eq!(
            spi.to_info().sampling_type,
            PufferfishType::SparseSparseV2(DEFAULT_CONTIG_SAMPLIG),
        );
    }

    #[test]
    fn validate_small_txome_no_wm() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseSparseIndex::from_sparse_index(pi, 100000, 1000000);
        assert!(spi.validate())
    }

    #[test]
    fn validate_small_txome_all_wm() {
        let p = to_abs_path(SMALL_TXOME_SPARSE_INDEX);
        let pi = SparseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseSparseIndex::from_sparse_index(pi, 0, 0);
        assert!(spi.validate())
    }
}
