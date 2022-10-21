use bincode::Result;
use clap::Parser;
use kmers::naive_impl::CanonicalKmer;
use rand::{rngs::StdRng, SeedableRng};
use serde::{Deserialize, Serialize};
use std::{fs::File, io::BufWriter, path::Path, time::Instant};

use pufferfish::benchmarking::{
    refseq::{RefSeq, RefSeqKmers},
    KmerSampler, RngKmer,
};
use pufferfish::index::sparse_contig_index::{StreamingCache, WMWalkCache};
use pufferfish::index::{prelude::*, QueryCache};

use crate::StatefulQuery;

use super::{IndexType, SparseIndexDriver, SparseSparseIndexDriver};

const SEED: u64 = 290348;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
pub struct GenKmerArgs {
    #[clap(short, long)]
    input: String,

    #[clap(short, long)]
    output: String,

    #[clap(short, long)]
    n_kmers: usize,

    #[clap(short, long, default_value_t = 31)]
    k: usize,

    #[clap(short, long, default_value_t=SEED)]
    seed: u64,
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
pub struct BenchKmerArgs {
    #[clap(short, long)]
    input_kmers: String,

    #[clap(subcommand)]
    index: IndexType,

    #[clap(short, long, default_value_t = 31)]
    k: usize,

    #[clap(short, long)]
    n_kmers: Option<usize>,
}

#[derive(Serialize, Deserialize)]
pub struct Kmers {
    pub kmers: Vec<u64>,
    pub k: usize,
}

pub fn run_bench_kmers(args: &BenchKmerArgs) -> Result<()> {
    log::info!("Loading bincode kmer set");
    let f = File::open(&args.input_kmers)?;
    let kmers: Vec<u64> = bincode::deserialize_from(f)?;
    let kmers = if let Some(n) = args.n_kmers {
        let n = usize::min(n, kmers.len());
        kmers[..n].to_owned()
    } else {
        kmers
    };

    let kmers = Kmers { kmers, k: args.k };

    match &args.index {
        IndexType::Sparse { index_dir } => {
            let pi = SparseIndex::deserialize_from_cpp(index_dir)?;
            let mut bencher = SparseIndexDriver {
                index: &pi,
                cache: QueryCache::new(),
            };
            bench_kmers(kmers, &mut bencher);
        }
        IndexType::SparseSparse {
            index_dir,
            k2u_dir,
            ctable_dir,
        } => {
            let pi = WMSparseSparseIndex::deserialize_from_parts(index_dir, k2u_dir, ctable_dir)?;
            let mut bencher = SparseSparseIndexDriver {
                index: &pi,
                sc: StreamingCache::new(),
                wc: WMWalkCache::new(super::WALK_CACHE_SIZE),
            };
            bench_kmers(kmers, &mut bencher);
        }
        _ => unimplemented!(),
    }
    Ok(())
}

fn bench_kmers<T: StatefulQuery>(kmers: Kmers, index: &mut T) {
    log::info!("Running benchmark w {}", index.index_name());
    let n_kmers = kmers.kmers.len();
    log::info!("Decoding {} kmer queries", n_kmers);

    let start_time = Instant::now();
    let k = kmers.k;
    for km in kmers.kmers {
        let km = CanonicalKmer::from_u64(km, k as u8);
        let hits = index.grp_no_cache(&km).expect("Must be true positive kmer");
        let mrps = index.as_vec_no_cache(&hits);
        criterion::black_box(mrps);
    }
    let elapsed = start_time.elapsed().as_secs_f32();
    log::info!("Result: queried {} kmers in {} seconds", n_kmers, elapsed);
}

// Bug: serialize the kmer class instead of vec of u64s.
pub fn run_gen_kmers(args: &GenKmerArgs) -> Result<()> {
    log::info!("Loading refseq from {}", args.input);
    let n_kmers = args.n_kmers;
    let start_time = Instant::now();
    let refseq = load_refseq(&args.input);
    let kmers = RefSeqKmers::new(refseq, args.k);
    let rng = StdRng::seed_from_u64(args.seed);
    let mut rng_kmers = KmerSampler::new(rng, kmers);
    log::info!("Generating {} random kmers to {}", n_kmers, args.output);
    let u64_kmers: Vec<u64> = (0..n_kmers).map(|_| rng_kmers.next_kmer_u64()).collect();

    let elapsed = start_time.elapsed().as_secs_f32();
    log::info!("Generated {} kmers in {} seconds", args.n_kmers, elapsed);
    log::info!("Saving to disk...");

    let f = File::create(&args.output)?;
    let writer = BufWriter::new(f);

    bincode::serialize_into(writer, &u64_kmers)?;

    Ok(())
}

fn load_refseq<P: AsRef<Path>>(dir: P) -> RefSeq {
    let p = dir.as_ref();
    RefSeq::deserialize_from_cpp(
        p.join("refseq.bin"),
        p.join("reflengths.bin"),
        p.join("refAccumLengths.bin"),
    )
    .unwrap()
}
