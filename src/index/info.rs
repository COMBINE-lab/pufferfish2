use super::ContigSamplingStrategy;
use super::PufferfishBase;
use super::PufferfishType;
use super::{
    BaseIndex, DenseIndex, 
    SparseIndex, 
    WMSparseDenseIndex, WMSparseSparseIndex,
};
use log::debug;
use serde::{Deserialize, Serialize};
use serde_json;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/******************************************************************************/
// Info (stored as json)
/******************************************************************************/

// FIXME: Split info to a legacy LegacyInfo or V1Info that contains info-fields
// that are have no Options. Fields like sample_size and extension_size should
// wholly be encoded in 'sampling_type'
// and 'sampling_type' should be renamed to pufferfish_type
#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct Info {
    pub index_version: u64,
    pub reference_gfa: Vec<String>,
    pub sampling_type: PufferfishType,

    #[serde(rename = "k")]
    pub kmer_size: u64,
    pub num_kmers: usize,
    pub num_contigs: usize,
    #[serde(rename = "seq_length")]
    pub seq_len: usize,
    pub have_ref_seq: bool,  //FIXME: has_ref_seq?
    pub have_edge_vec: bool, //FIXME: has_edge_vec?

    #[serde(rename = "SeqHash")]
    pub seq_hash: String,
    #[serde(rename = "NameHash")]
    pub name_hash: String,
    #[serde(rename = "SeqHash512")]
    pub seq_hash_512: String,
    #[serde(rename = "NameHash512")]
    pub name_hash_512: String,

    #[serde(rename = "DecoySeqHash")]
    pub decoy_seq_hash: String,
    #[serde(rename = "DecoyNameHash")]
    pub decoy_name_hash: String,
    pub num_decoys: usize,

    // only present if this is a sparse index
    pub sample_size: Option<usize>,
    pub extension_size: Option<usize>,

    pub first_decoy_index: usize,
    pub keep_duplicates: bool,
}
#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct BaseInfo {
    pub index_version: u64,
    pub reference_gfa: Vec<String>,
    #[serde(rename = "k")]
    pub kmer_size: u64,
    pub num_kmers: usize,
    pub num_contigs: usize,
    #[serde(rename = "seq_length")]
    pub seq_len: usize,
    pub has_ref_seq: bool,
    pub has_edge_vec: bool,

    #[serde(rename = "SeqHash")]
    pub seq_hash: String,
    #[serde(rename = "NameHash")]
    pub name_hash: String,
    #[serde(rename = "SeqHash512")]
    pub seq_hash_512: String,
    #[serde(rename = "NameHash512")]
    pub name_hash_512: String,

    #[serde(rename = "DecoySeqHash")]
    pub decoy_seq_hash: String,
    #[serde(rename = "DecoyNameHash")]
    pub decoy_name_hash: String,
    pub num_decoys: usize,

    pub first_decoy_index: usize,
    pub keep_duplicates: bool,
}

impl Info {
    pub fn load<P: AsRef<Path>>(p: P) -> Self {
        let p = p.as_ref();
        debug!("Loading Info from: {:?}", &p);
        let f = File::open(p).unwrap();

        let r = BufReader::new(f);

        serde_json::from_reader(r).unwrap()
    }

    // pub fn unwrap_sampling_strategy(&self) -> ContigSamplingStrategy {
    //     match self.sampling_type {
    //         PufferfishType::SparseContigDensePos(s) => s,
    //         PufferfishType::SparseContigSparsePos(s) => s,
    //         _ => panic!(
    //             "Cannot unwrap sampling strategy from {:?} pufferfish",
    //             self.sampling_type
    //         ),
    //     }
    // }
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum K2UInfo {
    Dense,
    Sparse {
        sample_size: usize,
        extension_size: usize,
    },
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum CTableInfo {
    Dense,
    Sparse(ContigSamplingStrategy),
    SparseWM {
        sampled_wm_thresh: usize,
        nonsamp_wm_thresh: usize,
        sampling: ContigSamplingStrategy,
    },
}

impl K2UInfo {
    // TODO: wtf is this name?
    pub fn to_dir_name(&self) -> String {
        match self {
            Self::Dense => "k2u_dense".to_string(),
            Self::Sparse {
                sample_size,
                extension_size,
            } => format!("k2u_sparse_s={}_e={}", sample_size, extension_size),
        }
    }
}

impl CTableInfo {
    pub fn to_dir_name(&self) -> String {
        match self {
            Self::Dense => "ctable_dense".to_string(),
            Self::SparseWM {
                sampled_wm_thresh,
                nonsamp_wm_thresh,
                sampling,
            } => {
                let strat_str = Self::contig_sampling_strat_to_dir_name(sampling);
                format!(
                    "ctable_sparseWM_s={}_n={}_{}",
                    sampled_wm_thresh, nonsamp_wm_thresh, strat_str
                )
            }

            Self::Sparse(strat) => {
                let strat_str = Self::contig_sampling_strat_to_dir_name(strat);
                format!("ctable_sparse_{}", strat_str)
            }
        }
    }

    pub fn contig_sampling_strat_to_dir_name(strat: &ContigSamplingStrategy) -> String {
        match strat {
            ContigSamplingStrategy::Minimal => "min".to_string(),
            ContigSamplingStrategy::Strided { stride } => format!("stride_s={}", stride),
            ContigSamplingStrategy::StridedBV { stride } => format!("stridebv_s={}", stride),
            ContigSamplingStrategy::StridedOccThresh { stride, thresh } => {
                format!("stride-samp-pop_s={}_t={}", stride, thresh)
            }
            ContigSamplingStrategy::StridedBVPop { stride, thresh } => {
                format!("stridebv-pop_s={}_t={}", stride, thresh)
            } // _ => todo!(),
        }
    }
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct V2Info {
    pub base_info: BaseInfo,
    pub k2u_info: K2UInfo,
    pub ctable_info: CTableInfo,
}

impl V2Info {
    pub fn load<P: AsRef<Path>>(p: P) -> Self {
        let p = p.as_ref();
        debug!("Loading Info from: {:?}", &p);
        let f = File::open(p).unwrap();

        let r = BufReader::new(f);

        serde_json::from_reader(r).unwrap()
    }

    pub fn from_dense_index(index: &DenseIndex) -> Self {
        let base_info = BaseInfo::from_base_index(index.get_base_index());

        Self {
            base_info,
            k2u_info: K2UInfo::Dense,
            ctable_info: CTableInfo::Dense,
        }
    }

    pub fn from_sparse_index(index: &SparseIndex) -> Self {
        let base_info = BaseInfo::from_base_index(index.get_base_index());

        Self {
            base_info,
            k2u_info: K2UInfo::Sparse {
                sample_size: index.sample_size,
                extension_size: index.extension_size,
            },
            ctable_info: CTableInfo::Dense,
        }
    }

    pub fn from_wm_sparse_sparse_index(index: &WMSparseSparseIndex) -> Self {
        let base_info = BaseInfo::from_base_index(index.get_base_index());

        Self {
            base_info,
            k2u_info: K2UInfo::Sparse {
                sample_size: index.sample_size,
                extension_size: index.extension_size,
            },
            ctable_info: CTableInfo::SparseWM {
                sampled_wm_thresh: index.sampled_wm_thresh,
                nonsamp_wm_thresh: index.nonsamp_wm_thresh,
                sampling: index.contig_sampling_strategy,
            },
        }
    }

    pub fn from_wm_sparse_dense_index(index: &WMSparseDenseIndex) -> Self {
        let base_info = BaseInfo::from_base_index(index.get_base_index());

        Self {
            base_info,
            k2u_info: K2UInfo::Dense,
            ctable_info: CTableInfo::SparseWM {
                sampled_wm_thresh: index.sampled_wm_thresh,
                nonsamp_wm_thresh: index.nonsamp_wm_thresh,
                sampling: index.contig_sampling_strategy,
            },
        }
    }
}

impl BaseInfo {
    pub fn from_base_index(index: &BaseIndex) -> Self {
        Self {
            index_version: index.index_version,
            reference_gfa: index.reference_gfa.clone(),
            kmer_size: index.k,
            num_kmers: index.num_kmers,
            num_contigs: index.num_contigs,
            seq_len: index.seq_len,
            has_ref_seq: index.has_ref_seq(),
            has_edge_vec: index.have_edge_vec,
            seq_hash: index.seq_hash.clone(),
            name_hash: index.name_hash.clone(),
            seq_hash_512: index.seq_hash_512.clone(),
            name_hash_512: index.name_hash_512.clone(),
            decoy_seq_hash: index.decoy_seq_hash.clone(),
            decoy_name_hash: index.decoy_name_hash.clone(),
            num_decoys: index.num_decoys,
            first_decoy_index: index.first_decoy_index,
            keep_duplicates: index.keep_duplicates,
        }
    }
}

impl BaseInfo {
    pub fn load<P: AsRef<Path>>(p: P) -> Self {
        let p = p.as_ref();
        debug!("Loading Info from: {:?}", &p);
        let f = File::open(p).unwrap();

        let r = BufReader::new(f);

        serde_json::from_reader(r).unwrap()
    }
}

impl K2UInfo {
    pub fn load<P: AsRef<Path>>(p: P) -> Self {
        let p = p.as_ref();
        debug!("Loading K2UInfo from: {:?}", &p);
        let f = File::open(p).unwrap();

        let r = BufReader::new(f);

        serde_json::from_reader(r).unwrap()
    }
}

impl CTableInfo {
    pub fn load<P: AsRef<Path>>(p: P) -> Self {
        let p = p.as_ref();
        debug!("Loading CTableInfo from: {:?}", &p);
        let f = File::open(p).unwrap();

        let r = BufReader::new(f);

        serde_json::from_reader(r).unwrap()
    }
}
