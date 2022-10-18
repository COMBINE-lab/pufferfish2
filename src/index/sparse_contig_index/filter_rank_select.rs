//filter_rank_select.rs
use core::ops::BitAnd;

// Convenience functions for filter, then rank and select for patterns

pub trait FilterRank<T: BitAnd> {
    // Given a filter `f`, index `i`, and slice
    // rank(i, f, arr) finds # elems where
    //     (arr[j] & f) == (f & arr[i])
    // for j in [0, i)
    // ranks begin at 0
    fn filter_rank(&self, i: usize, filter: T) -> usize;
}

pub trait FilterSelect<T: BitAnd> {
    // Given a filter and a pattern and rank `r`, finds the index where for the r-th time
    //     (arr[j] & f) == (f & p)
    // for j in [0, i)
    // ranks begin at 0, returns None if no elem matches (f & p)
    fn filter_select(&self, rank: usize, filter: T, pattern: T) -> Option<usize>;
}

impl<T> FilterRank<T> for [T]
where
    T: BitAnd<Output = T> + Eq + Copy,
{
    fn filter_rank(&self, i: usize, filter: T) -> usize {
        let pattern = filter & self[i];
        let slice = &self[..i];
        slice
            .iter()
            .take(i)
            .filter(|&x| (*x & filter) == pattern)
            .count()
    }
}

impl<T> FilterSelect<T> for [T]
where
    T: BitAnd<Output = T> + Eq + Copy,
{
    fn filter_select(&self, rank: usize, filter: T, pattern: T) -> Option<usize> {
        let pattern = filter & pattern;
        let mut accum = rank + 1;
        for (i, elem) in self.iter().enumerate() {
            if (*elem & filter) == pattern {
                accum -= 1
            }
            if accum == 0 {
                return Some(i);
            }
        }
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_filter_rank() {
        let s = &[true, false, true, false, false, true, true];
        let f = true;
        assert_eq!(s.filter_rank(0, f), 0);
        assert_eq!(s.filter_rank(1, f), 0);
        assert_eq!(s.filter_rank(2, f), 1);
        assert_eq!(s.filter_rank(3, f), 1);
        assert_eq!(s.filter_rank(4, f), 2);
        assert_eq!(s.filter_rank(5, f), 2);
        assert_eq!(s.filter_rank(6, f), 3);

        let s = &[0b000, 0b101, 0b111, 0b101, 0b010];
        assert_eq!(s.filter_rank(0, 0b111), 0);
        assert_eq!(s.filter_rank(3, 0b101), 2);
        assert_eq!(s.filter_rank(3, 0b111), 1);
        assert_eq!(s.filter_rank(4, 0b111), 0);
        assert_eq!(s.filter_rank(4, 0b101), 1);
        assert_eq!(s.filter_rank(3, 0b010), 2);
    }

    #[test]
    fn test_filter_sel() {
        let s = &[0b000, 0b101, 0b111, 0b101, 0b010];
        assert_eq!(s.filter_select(0, 0b101, 0b111), Some(1));
        assert_eq!(s.filter_select(1, 0b101, 0b111), Some(2));
        assert_eq!(s.filter_select(2, 0b101, 0b111), Some(3));
        assert_eq!(s.filter_select(0, 0b111, 0b111), Some(2));

        assert_eq!(s.filter_select(1, 0b111, 0b101), Some(3));
        assert_eq!(s.filter_select(1, 0b101, 0b000), Some(4));
        assert_eq!(s.filter_select(1, 0b101, 0b010), Some(4));

        assert_eq!(s.filter_select(2, 0b010, 0b101), Some(3));
        assert_eq!(s.filter_select(1, 0b010, 0b001), Some(1));
        assert_eq!(s.filter_select(0, 0b010, 0b000), Some(0));

        assert_eq!(s.filter_select(1, 0b101, 0b011), None);
    }

    #[test]
    fn filter_rank_sel_inverse() {
        let s = &[0b000, 0b101, 0b111, 0b101, 0b010];
        let idx = 3;
        let f = 0b010;
        let r = s.filter_rank(idx, f);
        assert_eq!(idx, s.filter_select(r, f, s[idx]).unwrap());
        assert_eq!(r, s.filter_rank(idx, f));
    }
}
