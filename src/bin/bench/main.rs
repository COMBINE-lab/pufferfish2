use ::kmers::naive_impl::CanonicalKmer;
use clap::Parser;
use std::fs::File;

use pufferfish::index::sparse_contig_index::StreamingCache;
use pufferfish::index::sparse_contig_index::WMWalkCache;
use pufferfish::index::{prelude::*, CachedQuery, QueryCache};
use pufferfish::log::{setup_file_logging, DEFAULT_LOG_FILE};

use bincode::Result;

mod kmers;
mod reads;

pub use reads::{bench_reads, Readset};

const WALK_CACHE_SIZE: usize = 5000;

// const SEED: u64 = 290348;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
enum CLI {
    GenKmers(kmers::GenKmerArgs),
    Kmers(kmers::BenchKmerArgs),

    Readset(reads::BenchReadsArgs),
    PrepFastq(reads::PrepReadsArgs),
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
pub enum IndexType {
    SparseSparse {
        #[clap(short, long)]
        index_dir: String,

        #[clap(short, long)]
        k2u_dir: String,

        #[clap(short, long)]
        ctable_dir: String,
    },
    SparseDense {
        #[clap(short, long)]
        index_dir: String,

        #[clap(short, long)]
        k2u_dir: String,

        #[clap(short, long)]
        ctable_dir: String,
    },
    Dense {
        index_dir: String,
    },
    Sparse {
        index_dir: String,
    },
}

fn main() -> Result<()> {
    setup_file_logging(DEFAULT_LOG_FILE);
    let args = CLI::parse();
    log::info!("Running {:?}", args);

    match &args {
        CLI::PrepFastq(args) => reads::prep_readset(args),
        CLI::Readset(args) => {
            log::info!("Loading bincode readset");
            let f = File::open(&args.input)?;
            let readset = bincode::deserialize_from(f)?;
            reads::run_bench_reads(&readset, &args.index)
        }

        CLI::GenKmers(args) => kmers::run_gen_kmers(args),
        CLI::Kmers(args) => kmers::run_bench_kmers(args),
    }
}

// Driver classes that own caches and pointers to immutable indices.
pub trait StatefulQuery {
    type HitsT: DecodeHit + std::fmt::Debug;
    fn grp(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT>;

    // No streaming cache
    fn grp_no_cache(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT>;

    fn as_vec(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos>;

    // No streaming cache
    fn as_vec_no_cache(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos>;

    fn k(&self) -> usize;
    fn index_name(&self) -> &str {
        std::any::type_name::<Self>()
    }
}

struct SparseIndexDriver<'a> {
    index: &'a SparseIndex,
    cache: QueryCache,
}

impl<'a> StatefulQuery for SparseIndexDriver<'a> {
    type HitsT = <SparseIndex as PuffQuery<'a>>::HitsT;
    fn grp(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.get_ref_pos_with_cache(km, &mut self.cache)
    }

    fn grp_no_cache(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.get_ref_pos(km)
    }

    fn as_vec(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec()
    }

    fn as_vec_no_cache(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec()
    }

    fn k(&self) -> usize {
        self.index.k()
    }
}

struct DenseIndexDriver<'a> {
    index: &'a DenseIndex,
    cache: QueryCache,
}

impl<'a> StatefulQuery for DenseIndexDriver<'a> {
    type HitsT = <DenseIndex as PuffQuery<'a>>::HitsT;
    fn grp(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.get_ref_pos_with_cache(km, &mut self.cache)
    }

    fn grp_no_cache(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.get_ref_pos(km)
    }

    fn as_vec(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec()
    }

    fn as_vec_no_cache(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec()
    }

    fn k(&self) -> usize {
        self.index.k()
    }
}

struct SparseSparseIndexDriver<'a> {
    index: &'a WMSparseSparseIndex,
    sc: StreamingCache,
    wc: WMWalkCache<'a, WMSparseSparseIndex>,
}

impl<'a> StatefulQuery for SparseSparseIndexDriver<'a> {
    type HitsT = <WMSparseSparseIndex as PuffQuery<'a>>::HitsT;
    fn grp(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.grp_w_cache(km, &mut self.sc)
    }

    fn grp_no_cache(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.get_ref_pos(km)
    }

    fn as_vec(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec_w_cache_cache(&mut self.wc, &mut self.sc)
    }

    fn as_vec_no_cache(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec_w_cache(&mut self.wc)
    }

    fn k(&self) -> usize {
        self.index.k()
    }
}

struct SparseDenseIndexDriver<'a> {
    index: &'a WMSparseDenseIndex,
    sc: StreamingCache,
    wc: WMWalkCache<'a, WMSparseDenseIndex>,
}

impl<'a> StatefulQuery for SparseDenseIndexDriver<'a> {
    type HitsT = <WMSparseDenseIndex as PuffQuery<'a>>::HitsT;
    fn grp(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.grp_w_cache(km, &mut self.sc)
    }

    fn grp_no_cache(&mut self, km: &CanonicalKmer) -> Option<Self::HitsT> {
        self.index.get_ref_pos(km)
    }

    fn as_vec(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec_w_cache_cache(&mut self.wc, &mut self.sc)
    }

    fn as_vec_no_cache(&mut self, hits: &Self::HitsT) -> Vec<MappedRefPos> {
        hits.as_vec_w_cache(&mut self.wc)
    }

    fn k(&self) -> usize {
        self.index.k()
    }
}
