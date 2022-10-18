pub mod contig_samplers;
pub mod filter_rank_select;

mod occs;

mod wm_ref_walker;
mod wm_sparse_contig_table;
mod wm_sparse_dense_index;
mod wm_sparse_projected_hits;
mod wm_sparse_sparse_index;

//Re-exports
pub use wm_ref_walker::WalkCache as WMWalkCache;

pub use wm_sparse_contig_table::SparseContigTable as WMSparseContigTable;
pub use wm_sparse_dense_index::SparseDenseIndex as WMSparseDenseIndex;
pub use wm_sparse_sparse_index::SparseSparseIndex as WMSparseSparseIndex;

use crate::index::contig::ContigOcc;
use contig_samplers::ContigSamplingStrategy;

pub const DEFAULT_CONTIG_SAMPLIG: ContigSamplingStrategy = ContigSamplingStrategy::Minimal;
pub struct StreamingCache {
    pub contig_start: usize,
    pub contig_len: usize,
    pub contig_id: usize,
    pub contig_occs: Option<Vec<ContigOcc>>,
}

impl StreamingCache {
    pub fn new() -> Self {
        Self {
            contig_start: usize::MAX,
            contig_len: usize::MAX,
            contig_id: usize::MAX,
            contig_occs: None,
        }
    }
}

impl Default for StreamingCache {
    fn default() -> Self {
        Self::new()
    }
}
