use serde::{Deserialize, Serialize};
use simple_sds::bit_vector::BitVector;
use simple_sds::bits::bit_len;
// use simple_sds::int_vector::IntVector;
use crate::serde_ext;
use simple_sds::ops::{BitVec, Rank, Select, SelectZero};

// WM only for byte size
type Symbol = u8;
const MAX_SYMBOL_WIDTH: usize = 8;

const BIT_MASK: [u8; 8] = [1, 2, 4, 8, 16, 32, 64, 128];

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct WaveletMatrixSlice<'a> {
    // for small alphabet size, we can just store the rank offsets...
    wm: &'a WaveletMatrix,
    start: usize,
    len: usize, // number of symbols in the WM
}

impl<'a> WaveletMatrixSlice<'a> {
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn access(&self, i: usize) -> Symbol {
        assert!(i < self.len);
        self.wm.access(self.start + i)
    }

    pub fn rank(&self, a: Symbol, i: usize) -> usize {
        assert!(i < self.len);
        self.wm.rank(a, self.start + i) - self.wm.rank(a, self.start)
    }

    pub fn select(&self, a: Symbol, r: usize) -> Option<usize> {
        let r = self.wm.rank(a, self.start) + r;
        Some(self.wm.select(a, r)? - self.start)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq)]
pub struct WaveletMatrix {
    // The wavelet matrix: An efficient wavelet tree for large alphabets
    // Claude, Navarro, Ordonez
    // Information Systems, 2015
    // https://doi.org/10.1016/j.is.2014.06.002
    #[serde(with = "serde_ext")]
    levels: Vec<BitVector>,
    zeros: Vec<usize>,         // "Z_l" array for WM
    alpha_offsets: Vec<usize>, // "C" array for extended WM

    len: usize,
    n_levels: usize,
    alpha_size: usize, // Alphabet size
    supports_rank: bool,
    supports_select: bool,
}

#[allow(clippy::len_without_is_empty)]
impl WaveletMatrix {
    pub fn get_slice(&self, start: usize, end: usize) -> WaveletMatrixSlice {
        let len = end - start;
        WaveletMatrixSlice {
            start,
            len,
            wm: self,
        }
    }
    pub fn with_ranksel_support(seq: &[Symbol], alpha_size: usize) -> Self {
        let mut wm = Self::new(seq, alpha_size);
        wm.enable_rank();
        wm.enable_select();
        wm
    }

    pub fn with_rank_support(seq: &[Symbol], alpha_size: usize) -> Self {
        let mut wm = Self::new(seq, alpha_size);
        wm.enable_rank();
        wm
    }

    pub fn with_select_support(seq: &[Symbol], alpha_size: usize) -> Self {
        let mut wm = Self::new(seq, alpha_size);
        wm.enable_select();
        wm
    }

    pub fn new(seq: &[Symbol], alpha_size: usize) -> Self {
        assert!(alpha_size > 1);
        Self::check_seq(seq, alpha_size);
        // simple_sds::bits::bit_len returns length of binary representation of
        // a given integer which is exactly ceil(log_2( ... ))
        let n_levels = bit_len((alpha_size - 1) as u64);

        assert!(n_levels <= bit_len(Symbol::MAX as u64));

        let len = seq.len();

        // Insert levels
        assert!(n_levels <= MAX_SYMBOL_WIDTH);
        let mut levels = Vec::with_capacity(n_levels);
        let mut zeros = Vec::with_capacity(n_levels);
        let mut seq = seq.to_vec();

        for &a in &seq {
            assert!((a as usize) < alpha_size);
        }

        for l in 0..n_levels {
            let mut n_zeros = 0;
            let mut lvl = Vec::with_capacity(len);
            let hi_bit = n_levels - 1 - l;

            for a in &seq {
                let bit = (a & BIT_MASK[hi_bit]) > 0;
                lvl.push(bit);
                if !bit {
                    n_zeros += 1;
                }
            }

            let bv = BitVector::from_iter(lvl.clone().into_iter());
            levels.push(bv);
            zeros.push(n_zeros);

            seq.sort_by_other(&lvl);
        }

        // Build "C" array.
        // If C[a] == usize::MAX, then a \notin C
        let mut alpha_offsets = vec![usize::MAX; alpha_size];

        if !seq.is_empty() {
            let first = seq[0];
            alpha_offsets[first as usize] = 0;
            for i in 1..seq.len() {
                let curr_a = seq[i];
                let prev_a = seq[i - 1];

                if prev_a != curr_a {
                    alpha_offsets[curr_a as usize] = i;
                }
            }
        }

        Self {
            levels,
            zeros,
            alpha_offsets,

            len,
            n_levels,
            alpha_size,

            supports_rank: false,
            supports_select: false,
        }
    }

    pub fn supports_access(&self) -> bool {
        self.supports_rank
    }

    pub fn supports_rank(&self) -> bool {
        self.supports_rank
    }

    pub fn supports_select(&self) -> bool {
        self.supports_select
    }

    pub fn enable_access(&mut self) {
        self.enable_rank();
    }

    pub fn enable_rank(&mut self) {
        self.supports_rank = true;
        for lvl in self.levels.iter_mut() {
            lvl.enable_rank();
        }
    }

    pub fn enable_select(&mut self) {
        self.supports_select = true;
        for lvl in self.levels.iter_mut() {
            lvl.enable_select();
            lvl.enable_select_zero();
        }
    }

    pub fn check_seq(seq: &[Symbol], alpha_size: usize) {
        let sigma = (alpha_size - 1) as Symbol;
        for &x in seq {
            assert!(x <= sigma, "{:#b} not allowed for sigma={}", x, sigma);
        }
    }

    pub fn height(&self) -> usize {
        self.n_levels
    }

    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn has_symbol(&self, a: Symbol) -> bool {
        self.alpha_offsets[a as usize] != usize::MAX
    }

    pub fn access(&self, i: usize) -> Symbol {
        assert!(self.supports_access());
        let mut j = i;
        let mut alpha = 0;

        for l in 0..self.n_levels {
            let bv = &self.levels[l];
            let hi_bit = self.n_levels - l - 1;

            if bv.get(j) {
                j = self.zeros[l] + bv.rank(j);
                alpha |= BIT_MASK[hi_bit];
            } else {
                j = bv.rank_zero(j);
            }
        }
        alpha
    }

    pub fn rank(&self, a: Symbol, i: usize) -> usize {
        assert!(self.supports_rank(), "Rank not supported");
        // returns the number of "a" in S[0, i)
        // NB: index-at-0 is different from WM paper, but same convention as simple_sds::BitVector
        // NB: rank_zero semantics for simple_sdsl assume infinate "0" padding to the rhs
        // So rank_zero can return arbitrarily large ranks from arbitrarily high indices
        // For WM, rank is only supported up to self.len()

        assert!((a as usize) < self.alpha_size);
        assert!(i <= self.len());

        if !self.has_symbol(a) {
            return 0;
        }

        let mut j = i;
        for l in 0..self.n_levels {
            let bv = &self.levels[l];
            let hi_bit = self.n_levels - 1 - l;
            if (a & BIT_MASK[hi_bit]) > 0 {
                j = self.zeros[l] + bv.rank(j);
            } else {
                j = bv.rank_zero(j);
            }
        }

        j - self.alpha_offsets[a as usize]
    }

    pub fn select(&self, a: Symbol, rank: usize) -> Option<usize> {
        assert!(self.supports_select(), "Select not supported");
        // Occurences start at 0.
        // select(a, rank(a, i)) == i

        if !self.has_symbol(a) {
            return None;
        }

        let mut j = self.alpha_offsets[a as usize] + rank;
        for l in (0..self.n_levels).rev() {
            let bv = &self.levels[l];
            let lo_bit = self.n_levels - 1 - l;
            if (a & BIT_MASK[lo_bit]) > 0 {
                j = bv.select(j - self.zeros[l])?;
            } else {
                j = bv.select_zero(j)?;
            }
        }
        Some(j)
    }
}

// From https://github.com/mbhall88/psdm/blob/0c8c4be5e4a6d566193b688824197fac2d233108/src/lib.rs#L13-L41
// And: https://internals.rust-lang.org/t/feature-request-methods-for-sorting-reordering-with-indices/15568/14
trait ArgSort<T: Ord> {
    fn argsort(&self) -> Vec<usize>;
}

trait SortExt<T: Ord> {
    fn reorder_unchecked(&mut self, order: &[usize]);
    fn sort_by_other<O: Ord>(&mut self, other: &[O]);
}

impl<T: Ord> ArgSort<T> for [T] {
    fn argsort(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..self.len()).collect();
        indices.sort_by_key(|&i| &self[i]);
        indices
    }
}

impl<T: Ord + Clone> SortExt<T> for Vec<T> {
    fn reorder_unchecked(&mut self, order: &[usize]) {
        // TODO: implement the in-place algorithm

        let seq = self.clone();
        for i in 0..order.len() {
            self[i] = seq[order[i]].clone()
        }
    }

    fn sort_by_other<O: Ord>(&mut self, other: &[O]) {
        let order = other.argsort();
        self.reorder_unchecked(&order);
    }
}

#[cfg(test)]
mod wm_tests {
    use super::WaveletMatrix as WM;
    use super::*;

    fn to_bv(v: Vec<u8>) -> BitVector {
        // map vec of u8s of 0s and 1s to a bv
        // makes it easy to make bvs from vec declarations like vec![ 0, 1, 1, 0, 1,]
        let bv = v.iter().map(|x| {
            assert!(*x <= 1);
            *x > 0
        });
        let mut bv = BitVector::from_iter(bv);
        bv.enable_rank();
        bv.enable_select();
        bv.enable_select_zero();
        bv
    }

    #[test]
    fn wm_len() {
        let symbols = vec![1, 2, 3];
        let wm = WM::with_ranksel_support(&symbols, 8);
        assert_eq!(wm.height(), 3);
        assert_eq!(wm.len(), 3);
    }

    #[test]
    fn wm_height() {
        let symbols = vec![1, 2, 3];
        let wm = WM::with_ranksel_support(&symbols, 8);
        assert_eq!(wm.height(), 3);

        let wm = WM::with_ranksel_support(&symbols, 7);
        assert_eq!(wm.height(), 3);

        let wm = WM::with_ranksel_support(&symbols, 9);
        assert_eq!(wm.height(), 4);

        let wm = WM::with_ranksel_support(&symbols, (u8::MAX as usize) + 1);
        assert_eq!(wm.height(), 8);
    }

    #[should_panic]
    #[test]
    fn wm_alpha_size_too_big() {
        let empty = Vec::new();
        WM::with_ranksel_support(&empty, (u8::MAX as usize) + 2);
    }

    #[should_panic]
    #[test]
    fn wm_alpha_size_too_small() {
        let empty = Vec::new();
        WM::with_ranksel_support(&empty, 1);
    }

    #[should_panic]
    #[test]
    fn check_seq_panics() {
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        WM::check_seq(&seq, 7);
    }

    #[test]
    fn wm_lvls() {
        // See Figure 4 in Claude et al. Information Systems 2015
        //  4   7   6   5   3   2   1   0   2   1   4   1   7
        // 100 111 110 101 011 010 001 000 010 001 100 001 111
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        let ll0 = vec![1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1];
        let ll1 = vec![1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1];
        let ll2 = vec![1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1];

        let ll0 = to_bv(ll0);
        let ll1 = to_bv(ll1);
        let ll2 = to_bv(ll2);

        let lvls = WM::with_ranksel_support(&seq, 8).levels;

        assert_eq!(lvls[0], ll0);
        assert_eq!(lvls[1], ll1);
        assert_eq!(lvls[2], ll2);
    }

    #[test]
    fn build_fig4() {
        // See Figure 4 in Claude et al. Information Systems 2015
        //  4   7   6   5   3   2   1   0   2   1   4   1   7
        // 100 111 110 101 011 010 001 000 010 001 100 001 111
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        let wm = WM::with_ranksel_support(&seq, 8);

        // Check WM height and length
        assert_eq!(wm.height(), 3);
        assert_eq!(wm.len(), seq.len());

        // Check levels
        let ll0 = vec![1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1];
        let ll1 = vec![1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1];
        let ll2 = vec![1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1];

        let ll0 = to_bv(ll0);
        let ll1 = to_bv(ll1);
        let ll2 = to_bv(ll2);

        let lvls = &wm.levels;

        assert_eq!(lvls[0], ll0);
        assert_eq!(lvls[1], ll1);
        assert_eq!(lvls[2], ll2);

        // Check Zs and Cs;
        assert_eq!(wm.zeros, vec![7, 7, 6]);

        // Last level encodes the ordering...
        // 0 4 4 2 2 6 1 1 1 5 3 7 7
        let c = vec![0, 6, 3, 10, 1, 9, 5, 11];
        assert_eq!(wm.alpha_offsets, c);
    }

    #[test]
    fn access_fig4() {
        // See Figure 4 in Claude et al. Information Systems 2015
        //  4   7   6   5   3   2   1   0   2   1   4   1   7
        // 100 111 110 101 011 010 001 000 010 001 100 001 111
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        let wm = WM::with_ranksel_support(&seq, 8);
        for (i, &s) in seq.iter().enumerate() {
            assert_eq!(s, wm.access(i));
        }
    }

    #[test]
    fn rank_fig4() {
        // See Figure 4 in Claude et al. Information Systems 2015
        // indices:
        //     0   1   2   3   4   5   6   7   8   9  10  11  12
        // vals:
        //     4   7   6   5   3   2   1   0   2   1   4   1   7
        //    100 111 110 101 011 010 001 000 010 001 100 001 111
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        let wm = WM::with_ranksel_support(&seq, 8);

        for a in 0..8 {
            let mut count = 0;
            for (i, &s) in seq.iter().enumerate() {
                assert_eq!(wm.rank(a, i), count);
                if s == a {
                    count += 1;
                }
            }
            assert_eq!(wm.rank(a, wm.len()), count);
        }
    }

    #[test]
    fn select_fig4() {
        // See Figure 4 in Claude et al. Information Systems 2015
        // indices:
        //     0   1   2   3   4   5   6   7   8   9  10  11  12
        // vals:
        //     4   7   6   5   3   2   1   0   2   1   4   1   7
        //    100 111 110 101 011 010 001 000 010 001 100 001 111
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        let wm = WM::with_ranksel_support(&seq, 8);
        for (i, &a) in seq.iter().enumerate() {
            let r = wm.rank(a, i);
            assert_eq!(wm.select(a, r), Some(i));
        }

        for a in 0..8 {
            assert_eq!(wm.select(a, 1000), None)
        }
    }

    #[test]
    #[should_panic]
    fn rank_panics() {
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        let wm = WM::with_ranksel_support(&seq, 8);
        wm.rank(0, 14);
    }

    #[test]
    #[should_panic]
    fn new_panics_alpha_too_big() {
        let seq = vec![4, 7, 6, 5, 3, 2, 1, 0, 2, 1, 4, 1, 7];
        let _ = WM::with_ranksel_support(&seq, 257);
    }

    #[test]
    fn ranksel_conventions() {
        let mut bv = BitVector::from_iter(vec![true, false, true].into_iter());
        bv.enable_rank();
        assert_eq!(bv.rank(0), 0);
        assert_eq!(bv.rank(2), 1);
        assert_eq!(bv.rank(100), 2);

        let mut bv = BitVector::from_iter(vec![true, false, false].into_iter());
        bv.enable_rank();
        assert_eq!(bv.rank(0), 0);
        assert_eq!(bv.rank(2), 1);
        assert_eq!(bv.rank(3), 1);
        assert_eq!(bv.rank(4), 1);

        // simple-sds and sdsl-lite semeantics assume infinite rhs 0s
        assert_eq!(bv.rank_zero(100), 99);
    }

    #[test]
    fn wm_empty_seq() {
        let seq = Vec::new();
        let sigma = 10;
        let wm = WM::with_ranksel_support(&seq, sigma);
        assert_eq!(wm.len(), 0);
        for a in 0..sigma {
            assert_eq!(wm.rank(a as u8, 0), 0);
        }

        for a in 0..sigma {
            assert_eq!(wm.select(a as u8, 0), None);
        }
    }

    #[test]
    fn wm_examples() {
        let examples = [
            (vec![0, 1, 2, 3, 0], 4),
            (vec![0, 1, 1, 1, 0, 1, 1, 1, 0], 2),
            (vec![0, 1, 2, 3, 8, 6, 7, 7, 1, 8, 9, 0], 10),
            (vec![0, 1, 2, 3, 8, 6, 7, 7, 1, 8, 9, 0], 256),
        ];

        for (seq, sigma) in examples {
            println!("{:?} {}", &seq, sigma);
            test_seq(&seq, sigma);
        }
    }

    // Larger tests
    fn test_seq(seq: &[Symbol], sigma: usize) {
        let wm = WM::with_ranksel_support(seq, sigma);
        assert_eq!(wm.height(), bit_len((sigma - 1) as u64));

        let sigma = sigma as Symbol;

        // access;
        let acc: Vec<Symbol> = (0..wm.len()).map(|i| wm.access(i)).collect();
        assert_eq!(seq, acc);

        // rank;
        for a in 0..sigma {
            let mut count = 0;
            for (i, &s) in seq.iter().enumerate() {
                assert_eq!(wm.rank(a, i), count, "(a:{}, i: {})", a, i);
                if s == a {
                    count += 1;
                }
            }
            assert_eq!(wm.rank(a, wm.len()), count);
        }

        // select;
        for (i, &a) in seq.iter().enumerate() {
            let r = wm.rank(a, i);
            assert_eq!(wm.select(a, r), Some(i));
        }

        for a in 0..sigma {
            assert_eq!(wm.select(a, wm.len()), None)
        }
    }
}

#[cfg(test)]
mod sort_tests {
    use super::*;
    #[test]
    fn argsort() {
        assert_eq!(vec!["c", "d", "b", "a"].argsort(), vec![3, 2, 0, 1]);
        let ll0 = vec![1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1];
        assert_eq!(
            ll0.argsort(),
            vec![4, 5, 6, 7, 8, 9, 11, 0, 1, 2, 3, 10, 12]
        );
    }

    #[test]
    fn reorder() {
        let mut v = vec!["c", "d", "b", "a"];
        let order = vec![2, 3, 1, 0];
        v.reorder_unchecked(&order);
        assert_eq!(v, vec!["b", "a", "d", "c"]);

        let mut a = vec![0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1];
        let b = vec![1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1];
        let o = vec![4, 5, 6, 7, 8, 9, 11, 0, 1, 2, 3, 10, 12];
        a.reorder_unchecked(&o);
        assert_eq!(a, b);
    }

    #[test]
    fn sort_by_other() {
        let mut v = vec!["c", "d", "b", "a"];
        let other = vec!["d", "b", "c", "a"];

        v.sort_by_other(&other);
        assert_eq!(v, vec!["a", "d", "b", "c"]);
    }
}
