use rayon::prelude::*;
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
        let mut idx = key_to_store as usize & self.mask;

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
        let mut idx = key_to_find as usize & self.mask;

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

fn main() {
    const K: usize = 21; // The table reserves one key bit, so K must be <= 31.

    // Allocate 2^24 slots upfront (~268MB of contiguous lock-free memory)
    let counter = LockFreeKmerCounter::new(24);

    // Simulate a DNA sequence.
    let sequence = b"ATGCTAGCTAGCTAGCTAGCGATCGATCGATCGATCGATCGATCGATC".repeat(200_000_000);

    // Parallel windows include k-mers that cross what would otherwise be
    // chunk boundaries.
    sequence.par_windows(K).for_each(|window| {
        let mut kmer_bits = 0u64;
        let mut valid = true;

        for &base in window {
            if let Some(bits) = base_to_bits(base) {
                kmer_bits = (kmer_bits << 2) | bits;
            } else {
                valid = false;
                break;
            }
        }

        if valid {
            counter
                .increment(kmer_bits)
                .expect("failed to increment k-mer counter");
        }
    });

    println!("K-mer counting complete. Threads exited successfully with zero lock overhead.");
}

#[cfg(test)]
mod tests {
    use super::LockFreeKmerCounter;
    use std::sync::Arc;
    use std::thread;

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
}
