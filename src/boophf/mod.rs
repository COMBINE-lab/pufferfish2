// boophf.rs
// Provides:
// - Load only compatibility with BooPHF.hpp
// - Does *not* provide functionality to build BooPHF
// - See: https://github.com/COMBINE-lab/pufferfish/blob/develop/include/BooPHF.hpp

use log::warn;
use serde::{Deserialize, Serialize};
use simple_sds::raw_vector::{AccessRaw, RawVector};
use std::cmp::Eq;
use std::collections::HashMap;
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::hash::Hash;
use std::io::Read;
use std::path::Path;

use crate::boophf::hash::{MultiHash, MultiHashResult, MultiHashState};
use crate::cpp::DeserializeFromCpp;
use crate::io::{read_u64_vec_with_len, ReadFrom};
use crate::serde_ext;

// Provide crate::booph::hash
pub mod hash;

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct BooPHF<T: Hash + Eq> {
    pub gamma: f64,
    //nb_levels: usize,
    last_bitset_rank: u64,
    pub n_elem: usize,
    pub levels: Vec<BoophfBitVec>,
    //final_hash_size: usize,
    final_hash: HashMap<T, u64>, //all hashes are u64
}

pub trait MPHF<T> {
    fn lookup(&self, item: &T) -> Option<u64>;
}

impl<T: Hash + Eq> PartialEq for BooPHF<T> {
    fn eq(&self, other: &Self) -> bool {
        warn!("Compare equality of BooPHFs at your own risk...");
        let mut res = self.last_bitset_rank == other.last_bitset_rank;
        res &= self.n_elem == other.n_elem;
        res &= self.levels == other.levels;
        res &= self.final_hash == other.final_hash;
        res
    }
}

impl<T: Hash + Eq> Eq for BooPHF<T> {}

impl<T: Hash + Eq + ReadFrom> ReadFrom for BooPHF<T> {
    fn read_from(f: &mut dyn Read) -> std::io::Result<Self> {
        let gamma = f64::read_from(f)?;
        let nb_levels = i32::read_from(f)? as usize;
        let last_bitset_rank = u64::read_from(f)?;
        let n_elem = u64::read_from(f)? as usize;

        let mut levels = Vec::with_capacity(nb_levels);

        for _ in 0..nb_levels {
            let bv = BoophfBitVec::read_from(f)?;
            levels.push(bv);
        }

        let final_hash_size = u64::read_from(f)? as usize;
        let mut final_hash = HashMap::new();

        for _ in 0..final_hash_size {
            let k = T::read_from(f)?;
            let v = u64::read_from(f)?;
            let result = final_hash.insert(k, v);
            assert!(result.is_none());
        }

        let mphf = BooPHF {
            gamma,
            // nb_levels,
            last_bitset_rank,
            n_elem,
            levels,
            // final_hash_size,
            final_hash,
        };

        Ok(mphf)
    }
}

impl<T: Hash + Eq + ReadFrom> DeserializeFromCpp for BooPHF<T> {
    fn deserialize_from_cpp<P: AsRef<Path>>(p: P) -> std::io::Result<Self> {
        let mut f = File::open(p)?;
        Self::read_from(&mut f)
    }
}

impl<T: Hash + Eq> MPHF<T> for BooPHF<T> {
    fn lookup(&self, item: &T) -> Option<u64> {
        let hash = self.lookup_in_levels(item);
        match hash {
            Some(_) => hash,
            None => self.lookup_in_final_hash(item),
        }
    }
}

impl<T: Hash + Eq> BooPHF<T> {
    #[inline]
    fn fast_range_64(word: u64, p: u64) -> u64 {
        // Map word into [0, p) *quickly* using
        // https://github.com/lemire/fastrange
        // (the implementation seems to be specified for u64 in BooPHF.hpp)
        let word = word as u128;
        let p = p as u128;
        let res = (word * p) >> 64;
        res as u64
    }

    fn lookup_in_levels(&self, item: &T) -> Option<u64> {
        // Look for item  in bloom filters
        let mut mh_state = MultiHashState::default();
        for (li, lvl) in self.levels.iter().enumerate() {
            // Not sure why the "statefulness" of multihash is *OUTSIDE* of the multihash
            let mh_res: MultiHashResult; // type annotation for *you* the reader.
            if li == 0 {
                mh_res = MultiHash::h0(item);
            } else if li == 1 {
                mh_res = MultiHash::h1(mh_state, item);
            } else {
                mh_res = MultiHash::next(mh_state);
            }
            mh_state = mh_res.state;
            // check if the bit according to multihash is set
            // Not sure why "hash_domain" is used in original impl...
            let hash = mh_res.hash;
            let pos = Self::fast_range_64(hash, lvl.len() as u64) as usize;
            let is_set = lvl.bit(pos);

            if is_set {
                // NB the ranks at each levels are offset by the size of levels
                // preceding it. So this "weird" func call is ok.
                let mphf_hash = lvl.rank(pos);
                return Some(mphf_hash as u64);
            }
        }

        None
    }

    fn lookup_in_final_hash(&self, item: &T) -> Option<u64> {
        let hash = self.final_hash.get(item);
        // Option.map<U,F>:= Option(F) -> Option(U)
        hash.map(|h| *h + self.last_bitset_rank)
    }

    pub fn collision_prob(&self) -> f64 {
        1. - self.miss_prob()
    }

    pub fn miss_prob(&self) -> f64 {
        self.levels.iter().map(|l| 1. - l.density()).product()
    }
}

/******************************************************************************/
// BitVec implementation for parity with:
// - level
// - bitvec
/******************************************************************************/
#[derive(Clone, Debug, Default, Serialize, Deserialize, Eq, PartialEq)]
pub struct BoophfBitVec {
    // number of bits
    // [removed] pub size: usize, // just use data.len()

    // number of machine words representing bits
    // just data.capacity() / sizeof(u64)
    // [removed] jpub n_words: usize,
    #[serde(with = "serde_ext")]
    pub data: RawVector, // Vec<u64> representation internally

    // rank of the first bit of each word
    // this can (and will) include an offset.
    // i.e. the rank of bit[0] can have rank > 0.
    pub ranks: Vec<usize>,
}

impl BoophfBitVec {
    const NB_BITS_PER_SAMPLE: usize = 512;

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn bit(&self, pos: usize) -> bool {
        // Match simple_sds AccessRaw Trait
        self.data.bit(pos)
    }

    pub fn count_ones(&self) -> usize {
        self.data.count_ones()
    }

    pub fn density(&self) -> f64 {
        (self.count_ones() as f64) / (self.len() as f64)
    }
}

// TODO instead Impl for simple_sds::ops::rank
trait Rank {
    fn rank(&self, pos: usize) -> usize;
}

impl Rank for BoophfBitVec {
    fn rank(&self, pos: usize) -> usize {
        let word_idx = pos / 64;
        let word_offset = pos % 64;
        let block = pos / Self::NB_BITS_PER_SAMPLE;

        let block_start_idx = block * Self::NB_BITS_PER_SAMPLE / 64;
        let mut r = self.ranks[block];
        for i in block_start_idx..word_idx {
            let word = self.data.word(i);
            r += word.count_ones() as usize;
        }
        let mask = (1 << word_offset) - 1;
        let offset = (self.data.word(word_idx) & mask).count_ones();
        r += offset as usize;
        r
    }
}

impl ReadFrom for BoophfBitVec {
    fn read_from(f: &mut dyn Read) -> std::io::Result<Self> {
        let n_bits = u64::read_from(f)? as usize;
        let n_words: usize = u64::read_from(f)? as usize;

        // NOTE
        // In CPP implementation of boophf, n_words is n_bits / 64 + 1
        // So there may be one extra serialized word if n_bits | 64.
        // RawVector::resize() truncates and fixes the issue
        let words = read_u64_vec_with_len(f, n_words)?;
        let data = RawVector::from_parts(n_bits, words);

        let ranks_size = u64::read_from(f)? as usize;
        let ranks = read_u64_vec_with_len(f, ranks_size)?;
        let ranks: Vec<usize> = ranks.iter().map(|x| *x as usize).collect();

        let bv = BoophfBitVec {
            // size,
            // n_words,
            data,
            ranks,
        };
        Ok(bv)
    }
}

impl<T: Hash + Eq> Display for BooPHF<T> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let omit = "[omitted]".to_string();
        f.debug_struct("BooPHF")
            .field("gamma", &self.gamma)
            // .field("nb_levels", &self.nb_levels)
            .field("last_bitset_rank", &self.last_bitset_rank)
            .field("n_elem", &self.n_elem)
            .field("levels", &omit)
            // .field("final_hash_size", &self.final_hash_size)
            .field("final_hash", &omit)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::fs::File;
    use std::path::Path;

    use simple_sds::bit_vector::BitVector;
    use simple_sds::ops::Rank as _;

    fn get_bbhash_small_example() -> BooPHF<u64> {
        // Load small example of bbhash produced by BooPHF.hpp
        // generated with first 10 uint64_t [0, 10)
        // with modified implementation that sets _nb_levels = 2
        let workdir = env!("CARGO_MANIFEST_DIR");
        let workdir = Path::new(workdir);
        let example_fp = "test_data/bbhash_n=10.bin";
        let example_fp = Path::new(example_fp);
        let example_fp = workdir.join(example_fp);

        let mut f = File::open(example_fp).unwrap();

        BooPHF::<u64>::read_from(&mut f).unwrap()
    }

    #[test]
    fn const_bits_per_block() {
        assert_eq!(BoophfBitVec::NB_BITS_PER_SAMPLE, 512);
    }

    #[test]
    fn parity_from_cpp_serialized() {
        let mphf = get_bbhash_small_example();
        assert_eq!(mphf.n_elem, 10);
        assert_eq!(mphf.final_hash.len(), 2);
        assert_eq!(mphf.levels.len(), 2);
    }

    #[test]
    fn test_level0_bits() {
        // ensure that the words are packed *exactly* the same as a
        // A verified C++ implementation
        let mphf = get_bbhash_small_example();
        let level0_word0 = mphf.levels[0].data.word(0);

        // Verified uint64_t representation of first packed word in BooPHF.hpp
        let cpp_level0_word0 = 2312599096050843650;
        assert_eq!(level0_word0, cpp_level0_word0);
    }

    #[test]
    fn test_levels_ranks() {
        // Make sure that BooPHF implementation matches
        // the ranks of sds implementation
        let mphf = get_bbhash_small_example();

        // boophf levels are offset by the rank of the last elem of
        // the preceding level
        let mut boophf_offset = 0;
        for bv in mphf.levels {
            // build an sds rank vector and compare
            let mut sds_bv = BitVector::from(bv.data.clone());
            sds_bv.enable_rank();

            for bit_pos in 0..bv.len() {
                let sds_rank = sds_bv.rank(bit_pos);
                let bb_rank = bv.rank(bit_pos);
                assert_eq!(sds_rank + boophf_offset, bb_rank);
                if bit_pos == (bv.len() - 1) {
                    boophf_offset += bb_rank;
                }
            }
        }
    }

    #[test]
    fn test_lookup0() {
        let mphf = get_bbhash_small_example();
        assert_eq!(mphf.lookup(&0), Some(2));
    }

    #[test]
    fn test_lookups() {
        let mphf = get_bbhash_small_example();
        let hashes: Vec<u64> = vec![2, 0, 8, 3, 5, 4, 1, 7, 6, 9, 7];
        for (item, checked) in hashes.iter().enumerate() {
            let item = item as u64;
            let hash = mphf.lookup(&item).unwrap();
            assert_eq!(hash, *checked);
        }
    }

    #[test]
    fn test_lookup_does_not_exist() {
        let mphf = get_bbhash_small_example();
        assert_eq!(mphf.lookup(&10), Some(7));
        assert_eq!(mphf.lookup(&11), Some(0));
        assert_eq!(mphf.lookup(&12), Some(0));
        assert_eq!(mphf.lookup(&13), None);

        for i in 13u64..20u64 {
            assert_eq!(mphf.lookup(&i), None);
        }
    }

    #[test]
    fn test_final_hash() {
        // we know that 2 and 9 map into the final level
        let mphf = get_bbhash_small_example();

        let hash9 = mphf.lookup_in_final_hash(&9);
        let hash2 = mphf.lookup_in_final_hash(&2);

        assert_eq!(hash9, Some(9));
        assert_eq!(hash2, Some(8));
    }
}
