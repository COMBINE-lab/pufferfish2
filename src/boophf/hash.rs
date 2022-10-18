// gen keys

//use crate::bitvector::BitVector;
use std::convert::TryInto;
use std::hash::{Hash, Hasher};
use std::num::Wrapping;

// Default seeds for simplehash and multihasher
const HASH_PAIR_SEEDS: (u64, u64) = (0xAAAAAAAA55555555, 0x33333333CCCCCCCC);

/******************************************************************************/
// SimpleHash
// -- For parity with SingleHashFunctor<uint64_t>
/******************************************************************************/

// SingleHashFunctor
pub struct SimpleHash {
    pub seed: u64,
    state: u64,
}

impl SimpleHash {
    fn new(seed: u64) -> Self {
        Self { seed, state: 0 }
    }

    pub fn hash_with_seed<T: Hash + ?Sized>(seed: u64, item: &T) -> u64 {
        let mut hasher = Self::new(seed);
        item.hash(&mut hasher);
        hasher.finish() // returns hash
    }

    fn hash64(key: u64, seed: u64) -> u64 {
        // allow overflow
        let mut hash = Wrapping(seed);
        let key = Wrapping(key);

        let init = (hash << 7) ^ (key * (hash >> 3)) ^ (!((hash << 11) + (key ^ (hash >> 5))));
        hash ^= init;
        hash = (!hash) + (hash << 21);
        hash = hash ^ (hash >> 24);
        hash = (hash + (hash << 3)) + (hash << 8);
        hash = hash ^ (hash >> 14);
        hash = (hash + (hash << 2)) + (hash << 4);
        hash = hash ^ (hash >> 28);
        hash = hash + (hash << 31);

        hash.0
    }
}

impl Hasher for SimpleHash {
    fn finish(&self) -> u64 {
        self.state
    }

    fn write(&mut self, bytes: &[u8]) {
        // TODO: a better hash function would consume byte-by-byte...
        if bytes.len() < std::mem::size_of::<u64>() {
            panic!("Do not use SimpleHash for < 8 byte data")
        }
        let (word, _rest) = bytes.split_at(std::mem::size_of::<u64>());
        //*input = rest; //don't need the rest of the bytes so no need to shift pointer
        let key = u64::from_ne_bytes(word.try_into().unwrap());
        //let key = u64::from_ne_bytes(bytes[..std::mem::size_of::<u64>()]);
        self.state = Self::hash64(key, self.seed);
    }
}

impl Default for SimpleHash {
    fn default() -> Self {
        Self {
            seed: HASH_PAIR_SEEDS.0,
            state: 0,
        }
    }
}
/******************************************************************************/
// MultiHash
// - For parity with XorshiftHashFunctors<u64, SingleHashFunctor>
// - typedef'd to MultiHasher_t in mphf
// - NB: This is kind of a stateful "chained" hash. Returning a tuple:
//       that is it's own state and the resulting hash value each time.
// - The implementation avoids mutable structs
/******************************************************************************/

pub struct MultiHashState(u64, u64);

impl Default for MultiHashState {
    fn default() -> Self {
        MultiHashState(HASH_PAIR_SEEDS.0, HASH_PAIR_SEEDS.1)
    }
}

pub struct MultiHashResult {
    pub state: MultiHashState,
    pub hash: u64,
}

// This probably should just be a module or keep it's own state
// but... for now have it own phantomdata
pub struct MultiHash<T> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: Hash + Sized> MultiHash<T> {
    // pub fn new() -> Self {
    //     Self {}
    // }

    pub fn h0(item: &T) -> MultiHashResult {
        let hash = SimpleHash::hash_with_seed(HASH_PAIR_SEEDS.0, &item);
        let state = MultiHashState(hash, HASH_PAIR_SEEDS.1);
        MultiHashResult { state, hash }
    }

    pub fn h1(state: MultiHashState, item: &T) -> MultiHashResult {
        let hash = SimpleHash::hash_with_seed(HASH_PAIR_SEEDS.1, &item);
        let state = MultiHashState(state.0, hash);
        MultiHashResult { state, hash }
    }
}
impl MultiHash<()> {
    pub fn next(state: MultiHashState) -> MultiHashResult {
        let s1 = std::num::Wrapping(state.0);
        let s0 = std::num::Wrapping(state.1);

        let s1 = s1 ^ (s1 << 23);
        let s1 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
        let hash = s1 + s0;
        let hash = hash.0; //unwrap after allowing for overflow
        let state = MultiHashState(s0.0, s1.0);

        MultiHashResult { state, hash }
    }
}

#[cfg(test)]
mod simple_hash_tests {
    use super::*;
    #[test]
    #[should_panic]
    fn panic_less_than_8_bytes() {
        SimpleHash::hash_with_seed(0, &0u32);
    }

    #[test]
    fn seeds() {
        assert_eq!(HASH_PAIR_SEEDS.0, SimpleHash::default().seed);
    }

    #[test]
    fn zero() {
        let state = SimpleHash::default();
        let hash_0 = 0x6e1bccdb7aa2bc25;
        let hash = SimpleHash::hash_with_seed(state.seed, &0u64);
        assert_eq!(hash_0, hash);
    }

    #[test]
    fn first10() {
        let state = SimpleHash::default();
        let true_hashes = vec![
            0x6e1bccdb7aa2bc25,
            0x54676a7b01425b7,
            0x5c9be323e5ad1be1,
            0x9567829f5e948f83,
            0xcf71e329165c79b5,
            0x9f1219f1bcd9d206,
            0x6bd828b35dba940e,
            0xf55b08c3340017c3,
            0xd178ae94742fa575,
            0x5dc299d49318dc6b,
        ];
        for (key, hash) in true_hashes.into_iter().enumerate() {
            let out = SimpleHash::hash_with_seed(state.seed, &key);
            assert_eq!(hash, out);
        }
    }
}

#[cfg(test)]
mod multi_hash_tests {
    use super::*;

    /* Check against known values produced from CPP binaries */
    #[test]
    fn zero_h0() {
        let key = &0u64;
        // let mh = MultiHash::new();
        let checked = 7934160411570650149;
        let res = MultiHash::h0(key);
        assert_eq!(checked, res.hash);
    }

    #[test]
    fn zero_h0_h1() {
        let key = &0u64;
        // let mh = MultiHash::new();

        let true_hashes = vec![7934160411570650149, 4031181471818755726];

        let res = MultiHash::h0(key);
        assert_eq!(true_hashes[0], res.hash);

        let res = MultiHash::h1(res.state, key);
        assert_eq!(true_hashes[1], res.hash);
    }

    #[test]
    fn zero_next() {
        // first "next" hash for zero
        let key = &0u64;
        // let mh = MultiHash::new();
        let res = MultiHash::h0(key);
        let res = MultiHash::h1(res.state, key);
        let res = MultiHash::next(res.state);
        let hash = res.hash;
        let checked = 7802733314557663513;
        assert_eq!(checked, hash);
    }

    #[test]
    fn zero_five() {
        // h0, h1, and first three "next" hashes
        let key = &0u64;
        // let mh = MultiHash::new();

        let true_hashes = vec![
            7934160411570650149,
            4031181471818755726,
            7802733314557663513,
            5772550616205298107,
            3882642898705877381,
        ];

        let mut impl_hashes = Vec::new();
        let res = MultiHash::h0(key);
        impl_hashes.push(res.hash);
        let res = MultiHash::h1(res.state, key);
        impl_hashes.push(res.hash);

        let mut state = res.state;
        for _ in 0..3 {
            let res = MultiHash::next(state);
            impl_hashes.push(res.hash);
            state = res.state;
        }

        for (cpp, rs) in true_hashes.iter().zip(impl_hashes.iter()) {
            assert_eq!(cpp, rs);
        }
    }
}
