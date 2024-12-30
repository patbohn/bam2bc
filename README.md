# Levenshtein distance barcode demultiplexing

This tool is used to demultiplex reads containing sequence-based barcodes. Specifically, it was developed for direct RNA sequencing, where barcodes were introduced following a constant sequence at the 3' end of the RNA. 

It calculates the minimal Levenshtein distance between the last N bases of the read and given barcode sequences. 

## How to install

Similar to other Rust tools, you first need to clone the source code (this repository). 

```bash
git clone https://github.com/patbohn/bam2bc
cd bam2bc
```

Then, assuming rust (and cargo) is set up, compile it:

```bash
cargo build --release
```

The binary will be at `bam2bc/target/release/bam2bc`

## How to run

Usage:
```bash
bam2bc 0.1.0
Detects sequence barcodes via Levenshtein distances given an adapter sequence.

USAGE:
    bam2bc [FLAGS] [OPTIONS] --barcode-file <barcode-file> -i <input-file>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v, --verbose    Will print to stdout

OPTIONS:
    -n, --adapter-search-region <adapter-search-region>
            number of nucleotides to search the adapter within at the ends [default: 200]

    -a, --barcode-file <barcode-file>
            file containing barcode sequences, one per row, now empty lines, no annotations

    -c, --constant-part <constant-part>
            constant part of sequence that is located 5' of the variable barcode [default: CATGGTCATAGCTGTTTCCTG]

    -e, --ignore-end-for-qscore <ignore-end-for-qscore>
            number of nucleotides to ignore at the 3' end to calculate mean qscore of read (ignore accidental DNA part
            basecalls) [default: 90]
    -i <input-file>                                        Input bam file
    -s, --max-samples-per-move <max-samples-per-move>
            Ignore moves with samples above this threshold to calculate mean translocation speed (correcting for
            stalling/tails) [default: 200]
    -m, --num-mismatches <num-mismatches>
            number of tolerated mismatches to find the constant sequence [default: 10]

    -o, --output-file <output-file>                        Output statistics file [default: output]
```
