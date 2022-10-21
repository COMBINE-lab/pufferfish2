# pufferfish2

`Pufferfish2` is a reference based index for exact k-mer queries designed to be a successor to
[`Pufferfish`](https://github.com/COMBINE-lab/pufferfish).

`Pufferfish2` sparsifies and compresses a `pufferfish` index by sampling unitigs and corresponding unitig-occurrences stored in `pufferfish`'s "unitig table", `utab`.

Currently, `pufferfish2` is also minimal reimplementation of `pufferfish` providing load-only compatibility.

## Building and testing

Compile `pufferfish2` with:

    cargo build --release

Run unit and integration tests with:

    cargo test

## Usage: building an index

After building a `pufferfish` index, build a `pufferfish2` index with:

    usage: sparsify_ctg_table [OPTIONS] --index-dir <INDEX_DIR> --out-index-dir <OUT_INDEX_DIR> <COMMAND>

    Commands:   
        strided-bv-pop      

    Options:
        -i, --index-dir <INDEX_DIR>                  
        -o, --out-index-dir <OUT_INDEX_DIR>          
        -s, --sampled-wm-thresh <SAMPLED_WM_THRESH>  [default: 64]
        -n, --nonsamp-wm-thresh <NONSAMP_WM_THRESH>  [default: 64]
        -h, --help                                   Print help information
        -V, --version                                Print version information

With subcommand:
    strided-bv-pop --stride <STRIDE> --thresh <THRESH>
    
    Options:
    -s, --stride <STRIDE>  
    -t, --thresh <THRESH>  
    -h, --help             Print help information

That samples positions of $\approx 1/s$ unitig and unitig occurrences, and all "popular" unitigs that occur more than `<THRESH>` times.


## Usage: benchmarking -- querying true positive kmers

NOTE: we currently only support `pufferfish` and `pufferfish2` indices built with "sparse" kmer-to-unitig mappings.

The provided benchmarking binaries first prepares input data and saves inputs to disks for repeatability and to avoid I/O bounds.

### Generate kmers
First generate k-mers from an existing *sparse* `pufferfish` index by pointing the binary to the directory containing a `pufferfish` index.
    Usage: bench gen-kmers [OPTIONS] --input <INPUT> --output <OUTPUT> --n-kmers <N_KMERS>

    Options:
    -i, --input <INPUT>      
    -o, --output <OUTPUT>    
    -n, --n-kmers <N_KMERS>  
    -k, --k <K>              [default: 31]
    -s, --seed <SEED>        [default: 290348]
    -h, --help               Print help information
    -V, --version            Print version information

### Run benchmark
Provide `<INPUT_KMERS>` generated kmer set, then provide to `pufferfish` or `pufferfish2` variants via `<COMMAND>`.
New `pufferfish2` indices using dense and sparse kmer-to-unitig mappings should be provided with `sparse-dense` and `sparse-sparse` commands, respectively. Old `pufferfish` variants using dense and sparse kmer-to-unitig mappings should be provided with `dense` and `sparse` commands, respectively.

    Usage: bench kmers [OPTIONS] --input-kmers <INPUT_KMERS> <COMMAND>

    Commands:
    sparse-sparse  
    sparse-dense   
    dense          
    sparse         
    help           Print this message or the help of the given subcommand(s)

    Options:
    -i, --input-kmers <INPUT_KMERS>  
    -k, --k <K>                      [default: 31]
    -n, --n-kmers <N_KMERS>          
    -h, --help                       Print help information
    -V, --version                    Print version information

## Usage: benchmarking -- querying reads from a FASTQ file

NOTE: we currently only support `pufferfish` and `pufferfish2` indices built with "sparse" kmer-to-unitig mappings.

### Prepare reads
Provide `<INPUT>` prepared readset, then provide to `pufferfish` or `pufferfish2` variants via `<COMMAND>`.


First load and prepare `<MAX_RADS>` reads from an input `.fastq` or `.fastq.gz` file --- `prep-fastq` processes and stores 2-bit encoded reads to disk.

    Usage: bench prep-fastq [OPTIONS] --input <INPUT> --output <OUTPUT>

    Options:
    -i, --input <INPUT>          
    -m, --max-reads <MAX_READS>  
    -o, --output <OUTPUT>        
    -h, --help                   Print help information
    -V, --version                Print version information

### Run benchmark

Provide `<INPUT>` prepared readset, then provide to `pufferfish` or `pufferfish2` variants via `<COMMAND>`.
New `pufferfish2` indices using dense and sparse kmer-to-unitig mappings should be provided with `sparse-dense` and `sparse-sparse` commands, respectively. Old `pufferfish` variants using dense and sparse kmer-to-unitig mappings should be provided with `dense` and `sparse` commands, respectively.

    Usage: bench readset --input <INPUT> <COMMAND>

    Commands:
    sparse-sparse  
    sparse-dense   
    dense          
    sparse         
    help           Print this message or the help of the given subcommand(s)

    Options:
    -i, --input <INPUT>  
    -h, --help           Print help information
    -V, --version        Print version information

## Micro-benchmarks

Run included micro-benchmarks (via `criterion`) for:
- implemented wavelet matrices

With:

    cargo bench