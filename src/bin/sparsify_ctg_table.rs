use clap::Parser;
use log::info;
use pufferfish::index::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
struct Args {
    #[clap(short, long)]
    index_dir: String,

    #[clap(short, long)]
    out_index_dir: String,

    #[clap(short, long, default_value_t = 64)]
    sampled_wm_thresh: usize,

    #[clap(short, long, default_value_t = 64)]
    nonsamp_wm_thresh: usize,

    #[clap(subcommand)]
    strategy: ContigSamplingStrategy,
}

fn main() {
    pufferfish::log::setup_default_logging();

    let args = Args::parse();

    let pb = PathBuf::from(&args.index_dir).join("info.json");
    let pinfo = Info::load(&pb);

    match pinfo.sampling_type {
        PufferfishType::Sparse => {
            info!("Sparsifying SparseIndex");
            info!("* from {}", args.index_dir);
            info!("* to {}", args.out_index_dir);

            info!("Loading SparseIndex");
            let load_time = Instant::now();
            let pi = SparseIndex::deserialize_from_cpp(&args.index_dir).unwrap();
            let load_time = load_time.elapsed();
            info!(
                "Finished loading SparseIndex in {:.3}s",
                (load_time.as_millis() as f64) * 1e-3
            );

            info!("Sparsifying contig table...");
            let spi = WMSparseSparseIndex::from_sparse_index_w_strategy(
                pi,
                args.strategy,
                args.sampled_wm_thresh,
                args.nonsamp_wm_thresh,
            );
            spi.lazy_serialize_to(args.out_index_dir).unwrap();
        }
        PufferfishType::Dense => {
            info!("Sparsifying DenseIndex");
            info!("* from {}", args.index_dir);
            info!("* to {}", args.out_index_dir);

            info!("Loading DenseIndex");
            let load_time = Instant::now();
            let pi = DenseIndex::deserialize_from_cpp(&args.index_dir).unwrap();
            let load_time = load_time.elapsed();
            info!(
                "Finished loading DenseIndex in {:.3}s",
                (load_time.as_millis() as f64) * 1e-3
            );

            info!("Sparsifying contig table...");
            let spi = WMSparseDenseIndex::from_dense_index_w_strategy(
                pi,
                args.strategy,
                args.sampled_wm_thresh,
                args.nonsamp_wm_thresh,
            );
            spi.lazy_serialize_to(args.out_index_dir).unwrap();
        }
        _ => panic!(
            "Cannot sparsify contig table of type {:?}",
            pinfo.sampling_type
        ),
    };
}
