use std::fs::File;
use std::io::{Read, Result};
use std::path::{Path, PathBuf};

use kmers::naive_impl::{Base, CanonicalKmer, Kmer};
use serde::{Deserialize, Serialize};
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Access, Vector};

use crate::cpp::{DeserializeFromCpp, FromCereal};

pub mod base_index;
pub mod compact;
pub mod consts;
pub mod contig;
mod dense_index;
pub mod info;
pub mod prelude;
pub mod projected_hits;
pub mod refseq_iter;
pub mod sparse_contig_index;
mod sparse_index;
pub mod validate;
// pub mod validate_fasta;

#[cfg(test)]
mod tests;

// Re-exports
pub use base_index::BaseIndex;
pub use contig::ContigOrientation;
pub use dense_index::DenseIndex;
pub use info::Info;
pub use sparse_contig_index::{WMSparseContigTable, WMSparseDenseIndex, WMSparseSparseIndex};
pub use sparse_index::SparseIndex;

use compact::FromCompact;
// use info::Info;
use refseq_iter::RefIterator;

use self::sparse_contig_index::contig_samplers::ContigSamplingStrategy;
use crate::seq::SeqVectorSlice;

// Pufferfish trait aliasing
// pub trait Pufferfish: for<'a> PuffQuery<'a> + PufferfishBase {}

// impl<'a> Pufferfish<'a> for DenseIndex {}
// impl<'a> Pufferfish<'a> for SparseIndex {}
pub trait Pufferfish: for<'a> PuffQuery<'a> + PufferfishBase {}
impl Pufferfish for DenseIndex {}
impl Pufferfish for SparseIndex {}
impl Pufferfish for WMSparseDenseIndex {}
impl Pufferfish for WMSparseSparseIndex {}

//****************************************************************************/
// Structs
//****************************************************************************/
#[derive(Debug, Eq, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub enum MappedOrientation {
    Forward,
    Backward,
}

#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize)]
pub struct MappedRefPos {
    pub ref_id: usize,
    pub pos: usize,
    pub orientation: MappedOrientation,
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct PosSamplingConfig {
    sample_size: usize,
    extension_size: usize,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone, Copy)]
pub enum PufferfishType {
    #[serde(rename = "dense")]
    Dense,
    #[serde(rename = "sparse")]
    Sparse,
    DenseRS,
    SparseRS,
    // SparseRS(PosSamplingConfig),
    SparseContigDensePos(ContigSamplingStrategy),
    SparseContigSparsePos(ContigSamplingStrategy),
    // SparseContigSparsePos(ContigSamplingStrategy, PosSamplingConfig),
    SparseDenseV2(ContigSamplingStrategy),
    SparseSparseV2(ContigSamplingStrategy),
}

impl PufferfishType {
    pub fn unwrap_contig_sampling_strategy(&self) -> ContigSamplingStrategy {
        match self {
            Self::SparseContigDensePos(s) => *s,
            Self::SparseContigSparsePos(s) => *s,
            Self::SparseDenseV2(s) => *s,
            Self::SparseSparseV2(s) => *s,
            _ => panic!("Cannot unwrap contig sampling strategy for {:?}", self),
        }
    }
}

pub struct QueryCache {
    pub contig_start: usize,
    pub contig_len: usize,
    pub contig_id: usize,
}

//****************************************************************************/
// Traits
//****************************************************************************/
pub trait DecodeHit {
    type IterT: Iterator<Item = MappedRefPos>;

    fn decode_hit(&self, i: usize) -> MappedRefPos;

    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
    fn contig_len(&self) -> usize;
    fn contig_id(&self) -> usize;
    fn mapped_orientation(&self) -> MappedOrientation;

    // Iter is most useful for early exit, returning Vec<MRP> can allow faster
    // compiler optimizations if "decode_hit(..)" requires branching

    fn iter(&self) -> Self::IterT;
    fn as_vec(&self) -> Vec<MappedRefPos>;
    // fn into_iter(self) -> Self::IterT
    // where
    //     Self: Sized; // maybe faster? some implementers use clone(...)
}

pub trait PuffQuery<'a> {
    type HitsT: DecodeHit + std::fmt::Debug;
    fn get_ref_pos(&'a self, kmer: &CanonicalKmer) -> Option<Self::HitsT>;
}

pub trait CachedQuery<'a, C> {
    type HitsT: DecodeHit + std::fmt::Debug;
    fn get_ref_pos_with_cache(
        &'a self,
        kmer: &CanonicalKmer,
        query_cache: &mut C,
    ) -> Option<Self::HitsT>;
}

pub trait PufferfishBase {
    fn get_encoded_contig_occs(&self, contig_id: usize) -> &[u64];

    fn get_ref_name(&self, id: usize) -> Option<&String>;
    fn get_ref_names(&self) -> &[String];
    fn to_info(&self) -> Info;

    fn get_base_index(&self) -> &BaseIndex;

    // Delegated to base index via default implementations
    // default functions are free when get_base_index is impl'd
    fn check_query(&self, kmer: &Kmer) {
        self.get_base_index().check_query(kmer)
    }

    fn seq_len(&self) -> usize {
        self.get_base_index().seq_len()
    }

    fn ref_len(&self, ref_id: usize) -> usize {
        self.get_base_index().ref_len(ref_id)
    }

    fn num_contigs(&self) -> usize {
        self.get_base_index().num_contigs()
    }
    fn num_kmers(&self) -> usize {
        self.get_base_index().num_kmers()
    }

    fn num_refs(&self) -> usize {
        self.get_base_index().num_refs()
    }

    fn num_total_occs(&self) -> usize;
    fn num_ctg_occs(&self, contig_id: usize) -> usize;

    fn get_ref_accum_lens(&self) -> &[u64] {
        self.get_base_index().get_ref_accum_lens()
    }
    fn k(&self) -> usize {
        self.get_base_index().k()
    }

    // Useful functions
    fn contig_len(&self, contig_id: usize) -> usize {
        self.get_base_index().contig_len(contig_id)
    }

    fn contig_start_pos(&self, contig_id: usize) -> usize {
        self.get_base_index().contig_start_pos(contig_id)
    }
    fn contig_prefix(&self, contig_id: usize, o: ContigOrientation) -> Kmer {
        self.get_base_index().contig_prefix(contig_id, o)
    }

    // Get kmers from global concatenated reference sequenc
    fn has_ref_seq(&self) -> bool {
        self.get_base_index().has_ref_seq()
    }

    fn get_refseq(&self, ref_id: usize) -> SeqVectorSlice {
        self.get_base_index().get_refseq(ref_id)
    }

    fn get_refseq_kmer_u64(&self, ref_id: usize, pos: usize) -> u64 {
        self.get_refseq(ref_id).get_kmer_u64(pos, self.k() as u8)
    }

    fn get_refseq_nuc(&self, ref_id: usize, pos: usize) -> Base {
        self.get_refseq(ref_id).get_base(pos)
    }
    fn get_refseq_kmer(&self, ref_id: usize, pos: usize) -> Kmer {
        let km = self.get_refseq_kmer_u64(ref_id, pos);
        Kmer::from_u64(km, self.k() as u8)
    }

    fn get_contig_nuc(&self, ctg_id: usize, ctg_o: ContigOrientation, offset: usize) -> Base {
        self.get_base_index().get_contig_nuc(ctg_id, ctg_o, offset)
    }

    fn get_contig_fw_kmer_u64(&self, ctg_id: usize, offset: usize) -> Base {
        self.get_base_index().get_contig_fw_kmer_u64(ctg_id, offset)
    }

    fn iter_refs(&self) -> RefIterator<Self>
    // fn iter_refs<'a, 'b>(&'a self) -> RefIterator<'b, Self>
    where
        Self: Sized;
    // fn get_ref_pos(&self, kmer: &CanonicalKmer) -> Option<ProjectedHits>;
}

#[derive(Clone)]
pub struct DenseContigTable {
    ctable: Vec<u64>,
    contig_offsets: IntVector,
    ref_names: Vec<String>,
    _ref_exts: Vec<u32>,
}

impl MappedRefPos {
    pub fn new(ref_id: usize, pos: usize, orientation: MappedOrientation) -> Self {
        MappedRefPos {
            ref_id,
            pos,
            orientation,
        }
    }

    pub fn is_forward(&self) -> bool {
        self.orientation.is_forward()
    }

    pub fn is_backward(&self) -> bool {
        !self.is_forward()
    }
}

impl MappedOrientation {
    pub fn reverse(&self) -> Self {
        match self {
            MappedOrientation::Forward => MappedOrientation::Backward,
            _ => MappedOrientation::Forward,
        }
    }

    pub fn is_forward(self) -> bool {
        self == MappedOrientation::Forward
    }

    pub fn is_backward(self) -> bool {
        !self.is_forward()
    }
}

impl DeserializeFromCpp for DenseContigTable {
    fn deserialize_from_cpp<P: AsRef<Path>>(dpath: P) -> Result<Self> {
        let files = PufferfishFilePaths::new(dpath);
        let (ref_names, ref_exts, ctable) = Self::deserialize_from_cpp_helper(files.ctable)?;
        let contig_offsets = IntVector::from_compact_serialized(files.ctg_offsets)?;
        Ok(Self {
            ctable,
            contig_offsets,
            ref_names,
            _ref_exts: ref_exts,
        })
    }
}

impl DenseContigTable {
    pub fn get_encoded_contig_occs(&self, ctg_id: usize) -> &[u64] {
        let start_idx = self.contig_offsets.get(ctg_id) as usize;
        let end_idx = self.contig_offsets.get(ctg_id + 1) as usize;
        let r = start_idx..end_idx;
        return self.ctable.get(r).unwrap();
    }

    pub fn num_ctg_occs(&self, ctg_id: usize) -> usize {
        let s = self.contig_offsets.get(ctg_id) as usize;
        let e = self.contig_offsets.get(ctg_id + 1) as usize;
        e - s
    }

    pub fn num_total_occs(&self) -> usize {
        self.ctable.len()
    }

    pub fn num_contigs(&self) -> usize {
        self.contig_offsets.len() - 1
    }

    pub fn get_offsets(&self) -> &IntVector {
        &self.contig_offsets
    }

    fn get_ref_name(&self, id: usize) -> Option<&String> {
        self.ref_names.get(id)
    }

    fn get_ref_names(&self) -> &[String] {
        &self.ref_names[..]
    }

    /**************************************************************************/
    // Helpers
    /**************************************************************************/
    fn deserialize_from_cpp_helper<P: AsRef<Path>>(
        p: P,
    ) -> std::io::Result<(Vec<String>, Vec<u32>, Vec<u64>)> {
        // Helper fn to deserialize serialized bits from cpp impl
        // that stores info as a tuple of vectors
        let mut f = File::open(p)?;

        // Deserialize Cereal input archive for Vec of std::string
        let ref_names = Vec::read_from_cereal_archive(&mut f)?;

        // Deserialize Cereal input archive for Vec of uint32
        let ref_exts = Vec::read_from_cereal_archive(&mut f)?;

        // Deserialize Cereal input archive for Vec of uint64
        let ctable = Vec::read_from_cereal_archive(&mut f)?;

        let mut buf = Vec::new();
        let n_bytes = f.read_to_end(&mut buf).unwrap();
        assert_eq!(n_bytes, 0);

        Ok((ref_names, ref_exts, ctable))
    }
}

impl QueryCache {
    pub fn new() -> Self {
        Self {
            contig_start: usize::MAX,
            contig_len: usize::MAX,
            contig_id: usize::MAX,
        }
    }
}

impl Default for QueryCache {
    fn default() -> Self {
        Self::new()
    }
}

/******************************************************************************/
// PufferfishFiles
//
// Dense index outputs files:
//    complete_ref_lens.bin
//    ctable.bin
//    ctg_offsets.bin
//    duplicate_clusters.tsv
//    info.json
//    mphf.bin
//    pos.bin
//    rank.bin
//    refAccumLengths.bin
//    ref_indexing.log
//    reflengths.bin
//    refseq.bin
//    seq.bin
/******************************************************************************/
#[derive(Debug)]
pub struct PufferfishFilePaths {
    // TODO: Absolute paths to all files?

    // Parent directory
    pub prefix: PathBuf,

    // binary format files
    pub complete_ref_lens: PathBuf,
    pub ctable: PathBuf,
    pub ctg_offsets: PathBuf,
    pub mphf: PathBuf, // Serialized BooPHF
    pub pos: PathBuf,
    pub sample_pos: PathBuf,        // [sparse index] sampled positions
    pub presence: PathBuf,          // [sparse index]
    pub extension_lengths: PathBuf, // [sparse inndex]
    pub extension_bases: PathBuf,   // [sparse index]
    pub direction: PathBuf,         // [sparse index]
    pub canonical: PathBuf,         // [sparse index]
    pub rank: PathBuf,              // contig boundaries
    pub ref_accum_lens: PathBuf,
    pub ref_lens: PathBuf,
    pub ref_seq: PathBuf,
    pub seq: PathBuf,

    // Other files
    pub duplicate_clusters_tsv: PathBuf,
    pub info_json: PathBuf,
    pub ref_indexing_log: PathBuf,
}

impl PufferfishFilePaths {
    pub fn new<P: AsRef<Path>>(dir: P) -> Self {
        let dir = dir.as_ref();
        Self {
            prefix: PathBuf::from(dir),
            complete_ref_lens: dir.join(consts::fp::COMPLETE_REF_LENS),
            ctable: dir.join(consts::fp::CTABLE),
            ctg_offsets: dir.join(consts::fp::CTG_OFFSETS),
            mphf: dir.join(consts::fp::MPHF),
            pos: dir.join(consts::fp::POS),
            sample_pos: dir.join(consts::fp::SAMPLE_POS),
            presence: dir.join(consts::fp::PRESENCE),
            extension_lengths: dir.join(consts::fp::EXTENSION_LENGTHS),
            extension_bases: dir.join(consts::fp::EXTENSION_BASES),
            direction: dir.join(consts::fp::DIRECTION),
            canonical: dir.join(consts::fp::CANONICAL),
            rank: dir.join(consts::fp::RANK),
            ref_accum_lens: dir.join(consts::fp::REF_ACCUM_LENS),
            ref_lens: dir.join(consts::fp::REF_LENS),
            ref_seq: dir.join(consts::fp::REF_SEQ),
            seq: dir.join(consts::fp::SEQ),
            duplicate_clusters_tsv: dir.join(consts::fp::DUPLICATE_CLUSTERS_TSV),
            info_json: dir.join(consts::fp::INFO_JSON),
            ref_indexing_log: dir.join(consts::fp::REF_INDEXING_LOG),
        }
    }

    // pub fn from_dir<P: AsRef<Path>>(dname: P) -> Self {
    //     let dir = dname.as_ref();
    //     // directory check is not checkable
    //     if !dir.exists() {
    //         panic!("The directory {} does not exist", dir.to_str().unwrap());
    //     }

    //     if !dir.is_dir() {
    //         panic!(
    //             "The path {} is not a valid directory",
    //             dir.to_str().unwrap()
    //         )
    //     }

    //     Self::with_dir(dname)
    // }
}

#[cfg(test)]
mod refseq_helper_tests {
    use super::dense_index::DenseIndex;
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::test_utils::*;
    use kmers::naive_impl::{Kmer, C, T};
    #[test]
    fn refseq_helpers() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        assert_eq!(pi.get_refseq_nuc(0, 5), T);
        assert_eq!(pi.get_refseq_nuc(1, 5), C);

        assert_eq!(pi.get_refseq_kmer(0, 1), Kmer::from("GTGAT"));
        assert_eq!(pi.get_refseq_kmer(1, 1), Kmer::from("GTGAC"));
    }
}

#[cfg(test)]
mod contig_helper_tests {
    use super::dense_index::DenseIndex;
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::test_utils::*;

    #[test]
    fn contig_len() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let ctg_a_kmer = CanonicalKmer::from("AGTGA");
        let ctg_b_kmer = CanonicalKmer::from("GTGAT");
        let ctg_c_kmer = CanonicalKmer::from("TGATA");
        let ctg_d_kmer = CanonicalKmer::from("GTAGA");
        let ctg_e_kmer = CanonicalKmer::from("AGGTA");
        let ctg_x_kmer = CanonicalKmer::from("GTAGC");
        let ctg_y_kmer = CanonicalKmer::from("GTAGC");
        let a_id = pi.get_ref_pos(&ctg_a_kmer).unwrap().contig_id;
        let b_id = pi.get_ref_pos(&ctg_b_kmer).unwrap().contig_id;
        let c_id = pi.get_ref_pos(&ctg_c_kmer).unwrap().contig_id;
        let d_id = pi.get_ref_pos(&ctg_d_kmer).unwrap().contig_id;
        let e_id = pi.get_ref_pos(&ctg_e_kmer).unwrap().contig_id;
        let x_id = pi.get_ref_pos(&ctg_x_kmer).unwrap().contig_id;
        let y_id = pi.get_ref_pos(&ctg_y_kmer).unwrap().contig_id;

        assert_eq!(pi.contig_len(a_id), 5);
        assert_eq!(pi.contig_len(b_id), 8);
        assert_eq!(pi.contig_len(c_id), 9);
        assert_eq!(pi.contig_len(d_id), 8);
        assert_eq!(pi.contig_len(e_id), 5);
        assert_eq!(pi.contig_len(x_id), 9);
        assert_eq!(pi.contig_len(y_id), 9);
    }

    #[test]
    fn contig_prefix() {
        let p = to_abs_path(TINY_REFS_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        let ctg_a_pref = CanonicalKmer::from("AGTGA");
        // let ctg_b_pref = CanonicalKmer::from("GTGAT");
        // let ctg_c_pref = CanonicalKmer::from("TGATA");
        let ctg_d_pref = CanonicalKmer::from("GTAGA");
        let ctg_e_pref = CanonicalKmer::from("AGGTA");
        // let ctg_x_pref = CanonicalKmer::from("GTAGC");
        // let ctg_y_pref = CanonicalKmer::from("GTAGC");
        let a_hits = pi.get_ref_pos(&ctg_a_pref).unwrap();
        let d_hits = pi.get_ref_pos(&ctg_d_pref).unwrap();
        let e_hits = pi.get_ref_pos(&ctg_e_pref).unwrap();

        // NB: cantigs are not written down in canonical sequence
        // let a_ctg_canonical = a_hits.mapped_orientation.is_forward();
        // let d_ctg_canonical = d_hits.mapped_orientation.is_forward();
        // let e_ctg_canonical = e_hits.mapped_orientation.is_forward();

        // // A contig is not canonical
        // // D contig is canonical
        // // E contig is not canonical

        // Canonical sequences of contigs
        // A: AGTGA
        // B: GTGATGAT
        // C: TGATAGTAG
        // D: GTAGAGGT
        // E: AGGTA
        // X: GTGACTGAT
        // Y: GTAGCAGGT

        // A
        assert_eq!(
            pi.contig_prefix(a_hits.contig_id, ContigOrientation::Forward),
            Kmer::from("TCACT")
        );

        assert_eq!(
            pi.contig_prefix(a_hits.contig_id, ContigOrientation::Backward),
            Kmer::from("AGTGA")
        );

        // D
        assert_eq!(
            pi.contig_prefix(d_hits.contig_id, ContigOrientation::Forward),
            Kmer::from("GTAGA")
        );

        assert_eq!(
            pi.contig_prefix(d_hits.contig_id, ContigOrientation::Backward),
            Kmer::from("GAGGT").to_reverse_complement()
        );

        // E
        assert_eq!(
            pi.contig_prefix(e_hits.contig_id, ContigOrientation::Forward),
            Kmer::from("AGGTA").to_reverse_complement()
        );

        assert_eq!(
            pi.contig_prefix(e_hits.contig_id, ContigOrientation::Backward),
            Kmer::from("AGGTA")
        );
    }
}
