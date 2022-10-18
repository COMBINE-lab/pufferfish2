use criterion::black_box;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use pufferfish::wm::WaveletMatrix as WM;
use rand::{rngs::StdRng, Rng, RngCore, SeedableRng};

// Thanks to this amazing blog:
// https://www.worthe-it.co.za/blog/2021-06-19-rust-performance-optimization-tools.html

const SEED: u64 = 2022;
const SIGMA: usize = 16;
const SIZES: &[usize] = &[50, 100, 500, 2000, 4000, 16000];

fn rand_seq<T: Rng>(rng: &mut T, len: usize, sigma: usize) -> Vec<u8> {
    let mut seq = vec![0; len];
    rng.fill_bytes(&mut seq);
    seq.iter().map(|&a| a % (sigma as u8)).collect()
}

fn linear_once(c: &mut Criterion) {
    let n = 1_000_000;
    let sigma = 16;

    let mut rng = StdRng::seed_from_u64(SEED);
    let seq = rand_seq(&mut rng, n, sigma);

    c.bench_function("linear_once", |b| {
        b.iter(|| {
            let a = rng.gen_range(0..sigma) as u8;
            let i = rng.gen_range(0..n);
            black_box(dumb_rank(&seq, a, i));
            // dumb_rank(black_box(&seq), a, i);
        })
    });
}

fn linear_once_unblackboxed(c: &mut Criterion) {
    let n = 1_000_000;
    let sigma = 16;

    let mut rng = StdRng::seed_from_u64(SEED);
    let seq = rand_seq(&mut rng, n, sigma);

    c.bench_function("linear_once_unblackboxed", |b| {
        b.iter(|| {
            let a = rng.gen_range(0..sigma) as u8;
            let i = rng.gen_range(0..n);
            dumb_rank(&seq, a, i);
            // dumb_rank(black_box(&seq), a, i);
        })
    });
}

fn wm_once(c: &mut Criterion) {
    let n = 1_000_000;
    let sigma = 16;

    let mut rng = StdRng::seed_from_u64(SEED);
    let seq = rand_seq(&mut rng, n, sigma);
    let wm = WM::with_ranksel_support(&seq, sigma);

    c.bench_function("wm_once", |b| {
        b.iter(|| {
            let a = rng.gen_range(0..sigma) as u8;
            let i = rng.gen_range(0..n);
            black_box(wm.rank(a, i));
        })
    });
}

fn wm_rank(c: &mut Criterion) {
    let mut group = c.benchmark_group("wm_rank");
    let mut rng = StdRng::seed_from_u64(SEED);
    for &len in SIZES {
        group.throughput(Throughput::Elements(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, &len| {
            // setup
            let seq = rand_seq(&mut rng, len, SIGMA);
            let wm = WM::with_ranksel_support(&seq, SIGMA);

            //query
            b.iter(|| {
                let a = rng.gen_range(0..SIGMA) as u8;
                let i = rng.gen_range(0..len);
                black_box(wm.rank(a, i));
            });
        });
    }
}

fn wm_ranksel(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut group = c.benchmark_group("wm_ranksel");
    for &len in SIZES {
        group.throughput(Throughput::Elements(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, &len| {
            // setup
            let seq = rand_seq(&mut rng, len, SIGMA);
            let wm = WM::with_ranksel_support(&seq, SIGMA);

            //query
            b.iter(|| {
                let i = rng.gen_range(0..len);
                let a = seq[i];
                let r = black_box(wm.rank(a, i));
                let ii = black_box(wm.select(a, r).unwrap());
                assert_eq!(i, ii);
            });
        });
    }
}

fn linear_rank(c: &mut Criterion) {
    let mut group = c.benchmark_group("linear_rank");
    let mut rng = StdRng::seed_from_u64(SEED);
    for &len in SIZES {
        group.throughput(Throughput::Elements(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, &len| {
            // setup
            let mut seq = vec![0; len];
            rng.fill_bytes(&mut seq);
            seq = seq.iter().map(|&a| a % (SIGMA as u8)).collect();

            //query
            b.iter(|| {
                let a = rng.gen_range(0..SIGMA) as u8;
                let i = rng.gen_range(0..len);
                black_box(dumb_rank(&seq, a, i));
            });
        });
    }
}

fn linear_rank_unblackboxed(c: &mut Criterion) {
    // This benchmark shows why blackboxing is necessary
    // if black_box is not used, then rustc can compile away "dumb_rank"
    let mut group = c.benchmark_group("linear_rank_unblackboxed");
    let mut rng = StdRng::seed_from_u64(SEED);
    for &len in SIZES {
        group.throughput(Throughput::Elements(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, &len| {
            // setup
            let mut seq = vec![0; len];
            rng.fill_bytes(&mut seq);
            seq = seq.iter().map(|&a| a % (SIGMA as u8)).collect();

            //query
            b.iter(|| {
                let a = rng.gen_range(0..SIGMA) as u8;
                let i = rng.gen_range(0..len);
                dumb_rank(&seq, a, i);
            });
        });
    }
}

fn linear_ranksel(c: &mut Criterion) {
    let mut group = c.benchmark_group("linear_ranksel");
    for &len in SIZES {
        group.throughput(Throughput::Elements(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, &len| {
            // setup
            let mut rng = StdRng::seed_from_u64(SEED);
            let mut seq = vec![0; len];
            rng.fill_bytes(&mut seq);
            seq = seq.iter().map(|&a| a % (SIGMA as u8)).collect();

            //query
            b.iter(|| {
                let i = rng.gen_range(0..len);
                let a = seq[i];
                let r = black_box(dumb_rank(&seq, a, i));
                let ii = black_box(dumb_select(&seq, a, r).unwrap());
                assert_eq!(i, ii);
            });
        });
    }
}

fn dumb_rank(seq: &[u8], a: u8, i: usize) -> usize {
    seq.iter().take(i).filter(|&x| *x == a).count()
}

fn dumb_select(seq: &[u8], a: u8, r: usize) -> Option<usize> {
    let mut count: usize = 0;
    for (i, &x) in seq.iter().enumerate() {
        if a == x {
            count += 1;
        }
        if count == (r + 1) {
            return Some(i);
        }
    }
    None
}

criterion_group! {
    name = wavelet_matrix;
    config = Criterion::default();
    targets = wm_once, linear_once, linear_once_unblackboxed, wm_rank, linear_rank, linear_rank_unblackboxed, wm_ranksel, linear_ranksel
    // targets = wm_ranksel, linear_ranksel,
}

criterion_main!(wavelet_matrix);
