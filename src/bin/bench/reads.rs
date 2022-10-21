use pufferfish::index::sparse_contig_index::{StreamingCache, WMWalkCache};
use pufferfish::index::{prelude::*, QueryCache};
use pufferfish::seq::SeqVector;

use kmers::naive_impl::CanonicalKmer;

use clap::Parser;
use criterion::black_box;
use log::{debug, info};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use std::time::Instant;

use super::StatefulQuery;
use bincode::Result;
use flate2::read::GzDecoder;

use super::{
    DenseIndexDriver, IndexType, SparseDenseIndexDriver, SparseIndexDriver,
    SparseSparseIndexDriver, WALK_CACHE_SIZE,
};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
pub struct PrepReadsArgs {
    #[clap(short, long)]
    pub input: String,

    #[clap(short, long)]
    pub max_reads: Option<usize>,

    #[clap(short, long)]
    pub output: String,
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
pub struct BenchReadsArgs {
    #[clap(short, long)]
    pub input: String,

    // #[clap(short, long)]
    // decode_thresh: Option<usize>,
    #[clap(subcommand)]
    pub index: IndexType,
}

pub fn prep_readset(args: &PrepReadsArgs) -> Result<()> {
    let outfile = &args.output;
    let readset = Readset::from_fastq_w_limit(&args.input, args.max_reads)?;

    let f = File::create(outfile)?;
    let writer = BufWriter::new(f);

    bincode::serialize_into(writer, &readset)
}

pub fn run_bench_reads(readset: &Readset, index_t: &IndexType) -> Result<()> {
    match index_t {
        IndexType::Dense { index_dir } => {
            let pi = DenseIndex::deserialize_from_cpp(index_dir)?;
            let mut bencher = DenseIndexDriver {
                index: &pi,
                cache: QueryCache::new(),
            };
            bench_reads(readset, &mut bencher);
        }

        IndexType::Sparse { index_dir } => {
            let pi = SparseIndex::deserialize_from_cpp(index_dir)?;
            let mut bencher = SparseIndexDriver {
                index: &pi,
                cache: QueryCache::new(),
            };
            bench_reads(readset, &mut bencher);
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
                wc: WMWalkCache::new(WALK_CACHE_SIZE),
            };
            bench_reads(readset, &mut bencher);
        }

        IndexType::SparseDense {
            index_dir,
            k2u_dir,
            ctable_dir,
        } => {
            let pi = WMSparseDenseIndex::deserialize_from_parts(index_dir, k2u_dir, ctable_dir)?;
            let mut bencher = SparseDenseIndexDriver {
                index: &pi,
                sc: StreamingCache::new(),
                wc: WMWalkCache::new(WALK_CACHE_SIZE),
            };
            bench_reads(readset, &mut bencher);
        }
    }
    Ok(())
}

#[derive(Serialize, Deserialize)]
pub struct Readset {
    reads: Vec<SeqVector>,
}

impl Readset {
    pub fn n_reads(&self) -> usize {
        self.reads.len()
    }
    pub fn n_kmers(&self, k: usize) -> usize {
        let mut n = 0;
        for read in &self.reads {
            n += read.len() - k + 1;
        }
        n
    }

    pub fn from_fastq_bufread<R: BufRead + std::fmt::Debug>(
        reader: &mut R,
        limit: Option<usize>,
    ) -> Result<Self> {
        let mut n_skipped = 0;
        let mut reads = Vec::new();

        let lines = reader.lines().skip(1).step_by(4);

        let denom = 1_000_000;
        for (i, line) in lines.enumerate() {
            if limit.is_some() && ((i - n_skipped) >= limit.unwrap()) {
                continue;
            }

            if i % denom == 0 {
                debug!("{}", i);
            }

            let line = line?;
            if line.contains('N') {
                // debug!("skipping record {}: {}", i, line);
                n_skipped += 1;
                continue;
            }

            let sv = SeqVector::from(&line);
            assert_eq!(String::from(&sv), line);
            reads.push(sv);
        }
        info!("Skipped {} reads with bad bases", n_skipped);
        Ok(Self { reads })
    }

    pub fn from_fastq_w_limit<P: AsRef<Path>>(fp: P, limit: Option<usize>) -> Result<Self> {
        let fp = fp.as_ref();
        let file = File::open(fp)?;

        if fp.extension().unwrap() == "gz" {
            debug!("Reading fastq.gz file...");
            let reader = GzDecoder::new(file);
            let mut reader = BufReader::new(reader);
            Self::from_fastq_bufread(&mut reader, limit)
        } else {
            debug!("Reading plain-text fastq file...");
            let mut reader = BufReader::new(file);
            Self::from_fastq_bufread(&mut reader, limit)
        }
    }

    pub fn from_fastq<P: AsRef<Path>>(fp: P) -> Result<Self> {
        Self::from_fastq_w_limit(fp, None)
    }

    pub fn iter_reads(&self) -> impl Iterator<Item = &SeqVector> {
        self.reads.iter()
    }
}

pub fn bench_reads<T: StatefulQuery>(readset: &Readset, bencher: &mut T) {
    info!("Running benchmark w {}", bencher.index_name());
    let mut n_kmers = 0;
    let mut n_kmer_hits = 0;
    let mut n_occs = 0;

    let k = bencher.k();
    info!("Num reads: {}", readset.n_reads());
    info!("Num Kmers: {}", readset.n_kmers(k));
    let start_time = Instant::now();

    for (ri, read) in readset.iter_reads().enumerate() {
        if ri % 1000 == 0 {
            log::debug!("read: {}", ri)
        }
        for km in read.iter_kmers(k as u8) {
            let km = CanonicalKmer::from(km);
            let result = bencher.grp(&km);
            if let Some(hits) = &result {
                n_occs += hits.len();
                n_kmer_hits += 1;
                let mrps = bencher.as_vec(hits);
                black_box(mrps);
            }
            n_kmers += 1;
        }
    }

    let elapsed = start_time.elapsed();
    info!("Num k_mer hits: {}", n_kmer_hits);
    info!(
        "Proportion of k-mers hit: {}",
        (n_kmer_hits as f64) / (n_kmers as f64)
    );
    info!("Num occurences found: {}", n_occs);
    info!("Elapsed Nanos: {}", elapsed.as_nanos());
    info!("Elapsed Micros: {}", elapsed.as_micros());
    info!("Elapsed Millis: {}", elapsed.as_millis());
    info!("Elapsed Seconds: {}", elapsed.as_secs_f64());
}
