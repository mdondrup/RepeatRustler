use bio::io::fasta;
use rayon::prelude::*;
use std::env;
use std::error::Error;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

// A single bucket in the Jellyfish-style hash table.
// Occupies exactly 16 bytes.
struct Bucket {
    // The high bit marks an occupied slot, leaving 63 bits for the k-mer key.
    kmer_key: AtomicU64,
    count: AtomicU64,
}

pub struct LockFreeKmerCounter {
    table: Vec<Bucket>,
    mask: usize, // Capacity is a power of two, so masking replaces modulo.
}

impl LockFreeKmerCounter {
    fn index(&self, mut key: u64) -> usize {
        // Mix all bases into the index instead of clustering shared suffixes.
        key = (key ^ (key >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        key = (key ^ (key >> 27)).wrapping_mul(0x94d049bb133111eb);
        ((key ^ (key >> 31)) as usize) & self.mask
    }

    /// Initializes a table with 2^power_of_two_size elements.
    pub fn new(power_of_two_size: u32) -> Self {
        let capacity = 1usize
            .checked_shl(power_of_two_size)
            .expect("table size exponent must be smaller than the usize bit width");
        let table = (0..capacity)
            .map(|_| Bucket {
                kmer_key: AtomicU64::new(0),
                count: AtomicU64::new(0),
            })
            .collect();

        Self {
            table,
            mask: capacity - 1,
        }
    }

    /// Adds one occurrence of a key containing at most 63 bits.
    ///
    /// The high bit is reserved to distinguish occupied buckets from empty
    /// buckets, so a 32-base k-mer cannot be represented by this table.
    pub fn increment(&self, key: u64) -> Result<(), &'static str> {
        const OCCUPIED_BIT: u64 = 1 << 63;
        if key >= OCCUPIED_BIT {
            return Err("k-mer key must fit in 63 bits");
        }

        let key_to_store = key | OCCUPIED_BIT;
        let mut idx = self.index(key);

        for _ in 0..self.table.len() {
            let bucket = &self.table[idx];
            let current_key = bucket.kmer_key.load(Ordering::Relaxed);

            if current_key == 0 {
                match bucket.kmer_key.compare_exchange(
                    0,
                    key_to_store,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => {
                        bucket.count.fetch_add(1, Ordering::Relaxed);
                        return Ok(());
                    }
                    Err(actual_key) if actual_key == key_to_store => {
                        bucket.count.fetch_add(1, Ordering::Relaxed);
                        return Ok(());
                    }
                    Err(_) => {}
                }
            } else if current_key == key_to_store {
                bucket.count.fetch_add(1, Ordering::Relaxed);
                return Ok(());
            }

            idx = (idx + 1) & self.mask;
        }

        Err("k-mer counter table is full")
    }

    #[cfg(test)]
    fn count(&self, key: u64) -> u64 {
        const OCCUPIED_BIT: u64 = 1 << 63;
        let key_to_find = key | OCCUPIED_BIT;
        let mut idx = self.index(key);

        for _ in 0..self.table.len() {
            let bucket = &self.table[idx];
            let current_key = bucket.kmer_key.load(Ordering::Relaxed);
            if current_key == 0 {
                return 0;
            }
            if current_key == key_to_find {
                return bucket.count.load(Ordering::Relaxed);
            }
            idx = (idx + 1) & self.mask;
        }
        0
    }

    fn write_text<W: Write>(&self, mut writer: W, k: usize) -> io::Result<()> {
        const OCCUPIED_BIT: u64 = 1 << 63;
        let mut entries = self
            .table
            .iter()
            .filter_map(|bucket| {
                let stored_key = bucket.kmer_key.load(Ordering::Relaxed);
                (stored_key & OCCUPIED_BIT != 0).then(|| {
                    (
                        stored_key & !OCCUPIED_BIT,
                        bucket.count.load(Ordering::Relaxed),
                    )
                })
            })
            .collect::<Vec<_>>();

        // Fixed-length two-bit keys have the same order as their DNA strings.
        entries.par_sort_unstable_by_key(|entry| entry.0);

        writeln!(writer, "kmer\tcount")?;
        let mut sequence = [b'A'; 31];
        for (key, count) in entries {
            decode_kmer(key, &mut sequence[..k]);
            writer.write_all(&sequence[..k])?;
            writeln!(writer, "\t{count}")?;
        }
        writer.flush()
    }
}

// Map characters to 2-bit sequences
#[inline(always)]
fn base_to_bits(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

fn decode_kmer(mut key: u64, sequence: &mut [u8]) {
    for base in sequence.iter_mut().rev() {
        *base = match key & 0b11 {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!(),
        };
        key >>= 2;
    }
}

fn count_sequence(
    sequence: &[u8],
    k: usize,
    counter: &LockFreeKmerCounter,
) -> Result<(), &'static str> {
    if !(1..=31).contains(&k) {
        return Err("k-mer length must be between 1 and 31");
    }
    if sequence.len() < k {
        return Ok(());
    }

    const CHUNK_SIZE: usize = 64 * 1024;
    let mask = (1u64 << (2 * k)) - 1;
    let windows = sequence.len() - k + 1;
    sequence[..windows]
        .par_chunks(CHUNK_SIZE)
        .enumerate()
        .try_for_each(|(chunk_index, starts)| {
            let start = chunk_index * CHUNK_SIZE;
            // Read the trailing overlap, but count each window in just one chunk.
            let bases = &sequence[start..start + starts.len() + k - 1];
            let mut key = 0u64;
            let mut valid_bases = 0usize;
            for &base in bases {
                if let Some(bits) = base_to_bits(base) {
                    key = ((key << 2) | bits) & mask;
                    valid_bases = (valid_bases + 1).min(k);
                    if valid_bases == k {
                        counter.increment(key)?;
                    }
                } else {
                    key = 0;
                    valid_bases = 0;
                }
            }
            Ok(())
        })
}

fn count_fasta_file(
    path: impl AsRef<Path>,
    k: usize,
    counter: &LockFreeKmerCounter,
) -> Result<(), Box<dyn Error>> {
    let reader = fasta::Reader::from_file(path.as_ref())?;
    for record in reader.records() {
        count_sequence(record?.seq(), k, counter).map_err(io::Error::other)?;
    }
    Ok(())
}

fn parse_cli_args<I>(args: I) -> Result<(OsString, OsString, usize), io::Error>
where
    I: IntoIterator<Item = OsString>,
{
    const DEFAULT_K: usize = 21;
    const USAGE: &str = "usage: repeat-rustler [-k <k>] <input.fasta> <output.tsv>";

    let mut k = DEFAULT_K;
    let mut input_path = None;
    let mut output_path = None;
    let mut iter = args.into_iter();

    while let Some(arg) = iter.next() {
        let text = arg.to_string_lossy();
        match text.as_ref() {
            "-k" | "--k" => {
                let value = iter.next().ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, USAGE))?;
                k = value
                    .to_string_lossy()
                    .parse::<usize>()
                    .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "k-mer size must be an integer"))?;
            }
            s if s.starts_with("-k=") => {
                let value = &s[3..];
                k = value
                    .parse::<usize>()
                    .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "k-mer size must be an integer"))?;
            }
            "-h" | "--help" => return Err(io::Error::new(io::ErrorKind::InvalidInput, USAGE)),
            _ => {
                if input_path.is_none() {
                    input_path = Some(arg);
                } else if output_path.is_none() {
                    output_path = Some(arg);
                } else {
                    return Err(io::Error::new(io::ErrorKind::InvalidInput, USAGE));
                }
            }
        }
    }

    let input_path = input_path.ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, USAGE))?;
    let output_path = output_path.ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, USAGE))?;
    Ok((input_path, output_path, k))
}

fn main() -> Result<(), Box<dyn Error>> {
    let (input_path, output_path, k) = parse_cli_args(env::args_os().skip(1))?;

    // Allocate 2^24 slots upfront (~268MB of contiguous lock-free memory).
    let counter = LockFreeKmerCounter::new(31);
    // Process one FASTA record at a time; windows never cross record boundaries.
    count_fasta_file(input_path, k, &counter)?;
    let output = BufWriter::new(File::create(output_path)?);
    counter.write_text(output, k)?;

    println!("K-mer counting complete. Results written to the output file.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{base_to_bits, count_fasta_file, count_sequence, parse_cli_args, LockFreeKmerCounter};
    use std::collections::{BTreeMap, BTreeSet};
    use std::ffi::OsString;
    use std::fs;
    use std::io;
    use std::sync::atomic::Ordering;
    use std::sync::Arc;
    use std::thread;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn parses_k_mer_size_flag_before_and_after_positional_args() {
        let parsed_before = parse_cli_args([
            OsString::from("-k"),
            OsString::from("17"),
            OsString::from("input.fasta"),
            OsString::from("output.tsv"),
        ])
        .unwrap();
        assert_eq!(parsed_before.2, 17);

        let parsed_after = parse_cli_args([
            OsString::from("input.fasta"),
            OsString::from("output.tsv"),
            OsString::from("--k"),
            OsString::from("31"),
        ])
        .unwrap();
        assert_eq!(parsed_after.2, 31);
    }

    #[test]
    fn counts_keys_without_conflating_zero_with_an_empty_slot() {
        let counter = LockFreeKmerCounter::new(2);

        counter.increment(0).unwrap();
        counter.increment(0).unwrap();
        counter.increment(1).unwrap();

        assert_eq!(counter.count(0), 2);
        assert_eq!(counter.count(1), 1);
        assert_eq!(counter.count(2), 0);
    }

    #[test]
    fn writes_sorted_kmer_counts_as_tab_separated_text() {
        let counter = LockFreeKmerCounter::new(3);
        counter.increment(0b11).unwrap(); // AT
        counter.increment(0b01).unwrap(); // AC
        counter.increment(0b01).unwrap(); // AC

        let mut output = Vec::new();
        counter.write_text(&mut output, 2).unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "kmer\tcount\nAC\t2\nAT\t1\n"
        );
    }

    #[test]
    fn reports_invalid_keys_and_a_full_table() {
        let counter = LockFreeKmerCounter::new(1);

        assert!(counter.increment(1 << 63).is_err());
        counter.increment(0).unwrap();
        counter.increment(1).unwrap();
        assert!(counter.increment(2).is_err());
    }

    #[test]
    fn increments_concurrently() {
        let counter = Arc::new(LockFreeKmerCounter::new(4));
        let threads: Vec<_> = (0..4)
            .map(|_| {
                let counter = Arc::clone(&counter);
                thread::spawn(move || {
                    for _ in 0..1_000 {
                        counter.increment(42).unwrap();
                    }
                })
            })
            .collect();

        for thread in threads {
            thread.join().unwrap();
        }

        assert_eq!(counter.count(42), 4_000);
    }

    #[test]
    fn rolling_counts_match_naive_windows_across_chunk_boundaries() {
        let mut state = 42u64;
        let mut sequence: Vec<_> = (0..(2 * 64 * 1024 + 100))
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGTacgt"[state as usize & 7]
            })
            .collect();
        for position in [0, 65_520, 65_535, 65_536, 65_550, 131_071] {
            sequence[position] = b'N';
        }
        sequence[100] = b'-';

        for k in [1, 2, 21, 31] {
            let mut expected = BTreeMap::<u64, u64>::new();
            for window in sequence.windows(k) {
                let key = window.iter().try_fold(0, |key, &base| {
                    base_to_bits(base).map(|bits| (key << 2) | bits)
                });
                if let Some(key) = key {
                    *expected.entry(key).or_default() += 1;
                }
            }

            let counter = LockFreeKmerCounter::new(18);
            count_sequence(&sequence, k, &counter).unwrap();
            for (&key, &count) in &expected {
                assert_eq!(counter.count(key), count, "k={k}, key={key}");
            }
            let total: u64 = counter
                .table
                .iter()
                .map(|bucket| bucket.count.load(Ordering::Relaxed))
                .sum();
            assert_eq!(total, expected.values().sum::<u64>(), "k={k}");
        }
    }

    #[test]
    fn handles_short_sequences_and_rejects_unsupported_lengths() {
        let counter = LockFreeKmerCounter::new(2);
        count_sequence(b"", 21, &counter).unwrap();
        count_sequence(b"ACGT", 21, &counter).unwrap();
        assert!(count_sequence(b"ACGT", 0, &counter).is_err());
        assert!(count_sequence(b"ACGT", 32, &counter).is_err());
        assert_eq!(counter.count(0), 0);
    }

    #[test]
    fn distributes_keys_with_identical_suffixes() {
        let counter = LockFreeKmerCounter::new(12);
        let indices: BTreeSet<_> = (0..4096).map(|key| counter.index(key << 24)).collect();
        assert!(indices.len() > 2000);
        for key in 0..4096 {
            counter.increment(key << 24).unwrap();
        }
        for key in 0..4096 {
            assert_eq!(counter.count(key << 24), 1);
        }
    }

    #[test]
    fn propagates_output_flush_errors() {
        struct FailingFlush;
        impl io::Write for FailingFlush {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                Ok(bytes.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Err(io::Error::other("flush failed"))
            }
        }

        let counter = LockFreeKmerCounter::new(1);
        assert!(counter.write_text(FailingFlush, 21).is_err());
    }

    #[test]
    fn reads_fasta_records_without_joining_record_boundaries() {
        let unique_id = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "repeat-rustler-{}-{unique_id}.fasta",
            std::process::id()
        ));
        fs::write(&path, b">first\nAAAA\n>second\nAAAA\n").unwrap();

        let counter = LockFreeKmerCounter::new(2);
        count_fasta_file(&path, 4, &counter).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(counter.count(0), 2);
    }
}
