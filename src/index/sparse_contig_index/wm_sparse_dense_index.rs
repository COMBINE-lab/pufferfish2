use bincode;

use log::debug;
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Access, Rank};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

use super::contig_samplers::ContigSamplingStrategy;
use super::occs::EncodedOccs;
use super::wm_sparse_contig_table::{SparseContigTable, SparseContigTableBuilder};
use super::wm_sparse_projected_hits::SparseProjHits;
use super::{StreamingCache, DEFAULT_CONTIG_SAMPLIG};

use crate::boophf::MPHF;
use crate::index::base_index::BaseIndex;
use crate::index::info::Info;
use crate::index::refseq_iter::RefIterator;
use crate::index::{DenseIndex, MappedOrientation, PuffQuery, PufferfishBase, PufferfishType};
use crate::io::{DeserializeFrom, SerializeTo};
use crate::serde_ext::{AsSerialize, FromDeserialize};
use kmers::naive_impl::{CanonicalKmer, MatchType};

#[derive(Debug, Eq, PartialEq)]
pub struct SparseDenseIndex {
    pos: IntVector,
    base: BaseIndex,
    ctg_table: SparseContigTable,

    // TODO pub crate?
    pub contig_sampling_strategy: ContigSamplingStrategy,

    pub sampled_wm_thresh: usize, // Threshhold # of occs to get WM support
    pub nonsamp_wm_thresh: usize, // Threshhold # of occs to get WM support
}

impl PufferfishBase for SparseDenseIndex {
    fn num_total_occs(&self) -> usize {
        self.ctg_table.num_total_occs()
    }
    fn num_ctg_occs(&self, contig_id: usize) -> usize {
        self.ctg_table.num_ctg_occs(contig_id)
    }

    fn get_encoded_contig_occs(&self, _contig_id: usize) -> &[u64] {
        panic!("Not available for WMSparseDenseIndex")
    }
    fn get_base_index(&self) -> &BaseIndex {
        &self.base
    }

    fn to_info(&self) -> Info {
        let base = &self.base;
        Info {
            sampling_type: PufferfishType::SparseDenseV2(self.contig_sampling_strategy),
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
            sample_size: None, // Putting into enum would break c++ compatiblity)
            extension_size: None,
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
    // fn get_ref_name(&self) -> String {  }

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

impl<'a> PuffQuery<'a> for SparseDenseIndex {
    type HitsT = SparseProjHits<'a, Self>;
    fn get_ref_pos(&self, kmer: &CanonicalKmer) -> Option<SparseProjHits<Self>> {
        self.check_query(&kmer.get_fw_mer());
        // position of match on seq
        let (pos, mapped_orientation) = self.get_kmer_pos(kmer)?;
        let pos = pos as usize;
        let contig_id = self.base.bv.rank(pos) as usize;

        let offset = self.contig_start_pos(contig_id);
        let contig_len = self.contig_len(contig_id);

        let offset = pos - offset;
        let hits = self.ctg_table.get_encoded_contig_occs(contig_id);

        Some(SparseProjHits {
            index: self,
            kmer: kmer.get_fw_mer(),
            mapped_orientation,
            offset,
            hits,
            contig_len,
            contig_id,
        })
    }
}

#[allow(dead_code)]
impl SparseDenseIndex {
    pub fn grp_w_cache(
        &self,
        kmer: &CanonicalKmer,
        cache: &mut StreamingCache,
    ) -> Option<SparseProjHits<Self>> {
        let (pos, mapped_orientation) = self.get_kmer_pos(kmer)?;
        let pos = pos as usize;
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
            kmer: kmer.get_fw_mer(),
            mapped_orientation,
            offset,
            hits: occs,
            contig_len,
            contig_id,
        };
        Some(hits)
    }

    pub fn from_dense_index_w_strategy(
        pi: DenseIndex,
        strategy: ContigSamplingStrategy,
        sampled_wm_thresh: usize,
        nonsamp_wm_thresh: usize,
    ) -> Self {
        // Can only build if ref_seq is available
        assert!(pi.base.ref_seq.is_some());
        let ctg_table =
            SparseContigTableBuilder::new(&pi, strategy, sampled_wm_thresh, nonsamp_wm_thresh)
                .build();
        Self {
            pos: pi.pos,
            base: pi.base,
            ctg_table,
            contig_sampling_strategy: strategy,
            sampled_wm_thresh,
            nonsamp_wm_thresh,
        }
    }

    pub fn from_dense_index(
        pi: DenseIndex,
        pop_sampled_thresh: usize,
        pop_nonsamp_thresh: usize,
    ) -> Self {
        Self::from_dense_index_w_strategy(
            pi,
            DEFAULT_CONTIG_SAMPLIG,
            pop_sampled_thresh,
            pop_nonsamp_thresh,
        )
    }

    fn get_kmer_pos(&self, kmer: &CanonicalKmer) -> Option<(u64, MappedOrientation)> {
        // return the position of queried kmer
        // on seq, concatenated sequence of contigs
        self.check_query(&kmer.get_fw_mer());
        let canonical_query = kmer.get_canonical_kmer();

        let word: u64 = canonical_query.into_u64();
        let hash = self.base.mphf.lookup(&word);

        match hash {
            None => None,
            Some(i) => {
                let candidate_pos = self.pos.get(i as usize);
                let kmer_at_pos = self
                    .base
                    .seq
                    .get_kmer(candidate_pos as usize, self.base.k as u8);

                match kmer.get_kmer_equivalency(&kmer_at_pos) {
                    MatchType::IdentityMatch => Some((candidate_pos, MappedOrientation::Forward)),
                    MatchType::TwinMatch => Some((candidate_pos, MappedOrientation::Backward)),
                    _ => None,
                }
            }
        }
    }

    pub fn num_sampled_contigs(&self) -> usize {
        self.ctg_table.num_sampled()
    }

    pub fn num_not_sampled_contigs(&self) -> usize {
        self.ctg_table.num_not_sampled()
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

        let v2_info = V2Info::from_wm_sparse_dense_index(self);

        let dir = dir.as_ref();
        let fps = SparseDenseIndexPaths::new_from_parts(
            dir,
            dir.join(v2_info.k2u_info.to_dir_name()),
            dir.join(v2_info.ctable_info.to_dir_name()),
        );

        self.serialize_base_lazy(fps.clone())?;
        self.serialize_ctable_lazy(fps.clone())?;
        self.serialize_k2u_lazy(fps)?;
        Ok(())
    }

    fn serialize_k2u_lazy(&self, fps: SparseDenseIndexPaths) -> bincode::Result<()> {
        let v2_info = V2Info::from_wm_sparse_dense_index(self);
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
            assert!(fps.pos.is_file());
        } else {
            debug!("Serializing K2U data structures");
            assert!(!fps.pos.is_file());
            std::fs::create_dir_all(fps.k2u_dir)?;

            let f: File = File::create(fps.k2u_info)?;
            serde_json::to_writer_pretty(f, &v2_info.k2u_info).unwrap();

            let f = File::create(fps.pos)?;
            let w = BufWriter::new(f);
            bincode::serialize_into(w, &self.pos.as_serialize())?;
        }
        Ok(())
    }

    fn serialize_base_lazy(&self, fps: SparseDenseIndexPaths) -> bincode::Result<()> {
        let v2_info = V2Info::from_wm_sparse_dense_index(self).base_info;

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

    fn serialize_ctable_lazy(&self, fps: SparseDenseIndexPaths) -> bincode::Result<()> {
        let v2_info = V2Info::from_wm_sparse_dense_index(self);
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
        let fps = SparseDenseIndexPaths::new_from_parts(dir, k2u_dir, ctable_dir);
        Self::deserialize_from_paths(&fps)
    }

    fn deserialize_from_paths(fps: &SparseDenseIndexPaths) -> bincode::Result<Self> {
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
        if let K2UInfo::Dense = k2u_info {
            // do nothing
        } else {
            panic!("Cannot deserialize WMSparseSparseIndex with {:?}", k2u_info)
        };

        debug!("Loading sparse contig table");
        let ctg_table = SparseContigTable::deserialize_from(&fps.ctable_dir)?;
        // ctg_table.enable_wm_ranksel();

        debug!("Loading pos");
        let f = File::open(&fps.pos)?;
        let r = BufReader::new(f);
        let pos = bincode::deserialize_from(r)?;
        let pos = IntVector::from_deserialized(pos);

        debug!("Loading base index");
        let f = File::open(&fps.base)?;
        let r = BufReader::new(f);
        let mut base: BaseIndex = bincode::deserialize_from(r)?;

        debug!("Loading MPHF");
        let f = File::open(&fps.mphf)?;
        let r = BufReader::new(f);
        let mphf = bincode::deserialize_from(r)?;
        base.mphf = mphf;

        let index = Self {
            pos,
            base,
            ctg_table,

            contig_sampling_strategy,
            sampled_wm_thresh,
            nonsamp_wm_thresh,
        };
        Ok(index)
    }
}
use crate::index::info::{BaseInfo, CTableInfo, K2UInfo, V2Info};

impl<P: AsRef<Path>> SerializeTo<P> for SparseDenseIndex {
    fn serialize_to(&self, dir: P) -> bincode::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir(dir)?;
        let fps = SparseDenseIndexPaths::new(dir);

        let v2_info = V2Info::from_wm_sparse_dense_index(self);
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

        let f = File::create(fps.pos)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.pos.as_serialize())?;

        Ok(())
    }
}

impl<P: AsRef<Path>> DeserializeFrom<P> for SparseDenseIndex {
    fn deserialize_from(dir: P) -> bincode::Result<Self> {
        let dir = dir.as_ref();
        let fps = SparseDenseIndexPaths::new(dir);
        Self::deserialize_from_paths(&fps)
    }
}

#[derive(Clone, Debug)]
pub struct SparseDenseIndexPaths {
    dir: PathBuf,
    k2u_dir: PathBuf,
    ctable_dir: PathBuf,

    info: PathBuf,
    k2u_info: PathBuf,
    ctable_info: PathBuf,

    pos: PathBuf,
    mphf: PathBuf,
    base: PathBuf,
}

impl SparseDenseIndexPaths {
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

            pos: k2u_dir.join("pos.bin"),
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
    // use crate::index::*;
    use crate::io::{DeserializeFrom, SerializeTo};
    use crate::test_utils::*;

    #[test]
    fn ser_de_tiny() {
        let dir = temp_file_name("");
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi, 0, 0);

        spi.serialize_to(&dir).unwrap();
        let spi_loaded = SparseDenseIndex::deserialize_from(&dir).unwrap();

        std::fs::remove_dir_all(dir).unwrap();
        assert_eq!(spi, spi_loaded);
    }
    #[test]
    fn ser_de_small() {
        let dir = temp_file_name("");
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi, 0, 0);

        spi.serialize_to(&dir).unwrap();
        let spi_loaded = SparseDenseIndex::deserialize_from(&dir).unwrap();

        std::fs::remove_dir_all(dir).unwrap();
        assert_eq!(spi, spi_loaded);
    }
}

#[cfg(test)]
mod wm_tests {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::index::DecodeHit;
    use crate::test_utils::*;

    #[test]
    fn pufferfish_type() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi, 0, 0);
        assert_eq!(
            spi.to_info().sampling_type,
            PufferfishType::SparseDenseV2(DEFAULT_CONTIG_SAMPLIG),
        );
    }
    #[test]
    fn tiny_multirefs_sampled_get_ref_pos() {
        // panic!();
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi.clone(), 0, 0);
        let ctg_a_kmer = CanonicalKmer::from("AGTGA");
        let ctg_e_kmer = CanonicalKmer::from("AGGTA");

        assert_iter_eq(
            spi.get_ref_pos(&ctg_a_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_a_kmer).unwrap().as_vec(),
        );
        assert_iter_eq(
            spi.get_ref_pos(&ctg_e_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_e_kmer).unwrap().as_vec(),
        );
    }

    #[test]
    fn tiny_multirefs_sampled_wm_slice() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi, 0, 0);
        let ctg_a_kmer = CanonicalKmer::from("AGTGA");
        let ctg_e_kmer = CanonicalKmer::from("AGGTA");

        let sph = spi.get_ref_pos(&ctg_a_kmer).unwrap();
        let occ = sph.hits.unwrap_sampled_wm();

        assert_eq!(occ.wm_slice.len(), 2);
        assert_eq!(occ.encoded_occs.len(), 2);
        assert_eq!(occ.len(), 2);

        // println!("{:?}", sph.hits);

        // println!("{:#4b}", occ.access_succ_o(0));
        // println!("{:#4b}", occ.access_succ_o(1));
        // println!("{:?}", SampledContigOcc::from(occ.encoded_occs[0]));
        // println!("{:?}", SampledContigOcc::from(occ.encoded_occs[1]));

        // 0-th occ is backwards, succeeded by C (0b00)
        // 1-st occ is backwards, succeeded by T (0b11)
        assert_eq!(occ.access_succ_o(0), 0b010);
        assert_eq!(occ.access_succ_o(1), 0b110);

        let sph = spi.get_ref_pos(&ctg_e_kmer).unwrap();
        let occ = sph.hits;
        let occ = occ.unwrap_sampled_wm();

        assert_eq!(occ.wm_slice.len(), 2);
        assert_eq!(occ.encoded_occs.len(), 2);
        assert_eq!(occ.len(), 2);

        // The last contig should have no successors
        // So the value at index i for each occ should be 0b1000
        // The sentinel "stub" value "0b1000"
        assert_eq!(occ.access_succ_o(0), 0b1000);
        assert_eq!(occ.access_succ_o(1), 0b1000);
    }

    #[test]
    fn tiny_multirefs_nonsampled_get_ref_pos_pred_sampled() {
        // Only query kmers whose on unitigs that immediately succeed a sampled unitig
        // setting popular threshold guarantees that sampled contigs are sampled not in
        // wavelet matrices. (This is the easiest not-sampled case)
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi.clone(), 100, 0);
        let ctg_b_kmer = CanonicalKmer::from("GTGAT");
        // let ctg_c_kmer = CanonicalKmer::from("TGATA");
        // let ctg_d_kmer = CanonicalKmer::from("GTAGA");

        // let ctg_x_kmer = CanonicalKmer::from("GTAGC");
        // let ctg_y_kmer = CanonicalKmer::from("GTAGC");

        assert_eq!(spi.get_ref_pos(&ctg_b_kmer).unwrap().len(), 1);

        assert_eq!(
            spi.get_ref_pos(&ctg_b_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_b_kmer).unwrap().as_vec(),
        );
    }

    #[test]
    fn tiny_multirefs_nonsampled_get_ref_pos_no_wm() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi.clone(), 100, 100);
        let ctg_b_kmer = CanonicalKmer::from("GTGAT");
        let ctg_c_kmer = CanonicalKmer::from("TGATA");
        let ctg_d_kmer = CanonicalKmer::from("GTAGA");

        let ctg_x_kmer = CanonicalKmer::from("GTAGC");
        let ctg_y_kmer = CanonicalKmer::from("GTAGC");

        assert_eq!(spi.get_ref_pos(&ctg_b_kmer).unwrap().len(), 1);

        assert_eq!(
            spi.get_ref_pos(&ctg_b_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_b_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_c_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_c_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_d_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_d_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_x_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_x_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_y_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_y_kmer).unwrap().as_vec(),
        );
    }

    #[test]
    fn tiny_multirefs_nonsampled_get_ref_pos_all_wm() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi.clone(), 0, 0);
        let ctg_b_kmer = CanonicalKmer::from("GTGAT");
        let ctg_c_kmer = CanonicalKmer::from("TGATA");
        let ctg_d_kmer = CanonicalKmer::from("GTAGA");

        let ctg_x_kmer = CanonicalKmer::from("GTAGC");
        let ctg_y_kmer = CanonicalKmer::from("GTAGC");

        assert_eq!(spi.get_ref_pos(&ctg_b_kmer).unwrap().len(), 1);

        assert_eq!(
            spi.get_ref_pos(&ctg_b_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_b_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_c_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_c_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_d_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_d_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_x_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_x_kmer).unwrap().as_vec(),
        );

        assert_eq!(
            spi.get_ref_pos(&ctg_y_kmer).unwrap().as_vec(),
            pi.get_ref_pos(&ctg_y_kmer).unwrap().as_vec(),
        );
    }

    use crate::index::validate::*;
    #[test]
    fn validate_small_txome_no_wm() {
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi, 100000, 1000000);
        assert!(spi.validate())
    }

    #[test]
    fn validate_small_txome_all_wm() {
        let p = to_abs_path(SMALL_TXOME_DENSE_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let spi = SparseDenseIndex::from_dense_index(pi, 0, 0);
        assert!(spi.validate())
    }
}
