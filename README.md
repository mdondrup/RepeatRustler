# RepeatRustler

Count 21-mers from a FASTA file and write the counts as a tab-separated text
file:

```sh
cargo run --release -- path/to/sequences.fasta kmer-counts.tsv
cargo run --release -- -k 31 path/to/sequences.fasta kmer-counts.tsv
```

The program streams FASTA records and counts valid DNA k-mers within each
record. Windows do not cross record boundaries, and bases outside A/C/G/T
invalidate any k-mer that contains them. The output has `kmer` and `count`
columns and is sorted by k-mer. `-k` sets the k-mer size and defaults to 21.

## Performance and comparison with Jellyfish

The counter uses rolling two-bit encoding, processes overlapping chunks in
parallel, and mixes the entire encoded key before probing the hash table.
Sorting uses numeric keys (equivalent to lexicographic order for fixed-length
k-mers), and DNA strings are decoded only while writing the output.
The fixed table has 2^24 slots; counting fails explicitly if it fills.

Build before timing so compilation is not included:

```sh
cargo build --release
time target/release/repeat-rustler sequences.fasta kmer-counts.tsv
```

For a comparable Jellyfish workload, use `-m 21`, not `-m 31`, and do not
enable canonical counting (`-C`): RepeatRustler counts each strand separately.
`jellyfish count` writes a binary database, whereas RepeatRustler also sorts
and writes a text table. Include Jellyfish's dump and sorting steps when
comparing end-to-end text output. Writing RepeatRustler's output to
`/dev/null` still includes sorting and text formatting, but avoids file storage.

On the local yeast input, single runs of the pre-change and optimized release
executables writing to `/dev/null` took 104.31 s and 0.63 s wall time,
respectively. The optimized run writing a TSV took 0.96 s and matched the
existing TSV byte-for-byte. These are local measurements, not a general
performance guarantee or a matched Jellyfish benchmark.