// Load only compatibility with cpp compact
use simple_sds::bit_vector::BitVector;
use simple_sds::int_vector::IntVector;

use simple_sds::raw_vector::RawVector;
use std::fs::File;
use std::io::Result;
use std::path::Path;

use crate::io::ReadFrom;
use crate::seq::SeqVector;

#[inline]
#[allow(dead_code)]
pub fn get_bits_per_element<P: AsRef<Path>>(p: P) -> usize {
    let mut f = File::open(p).unwrap();
    let _static_flag = u64::read_from(&mut f).unwrap();
    let width = u64::read_from(&mut f).unwrap();
    width as usize
}

pub trait FromCompact {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self>
    where
        Self: Sized;
}

impl FromCompact for SeqVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let iv = IntVector::from_compact_serialized(p)?;
        Ok(Self::from(iv))
    }
}

impl FromCompact for RawVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let iv = IntVector::from_compact_serialized(p)?;
        Ok(RawVector::from(iv))
    }
}

impl FromCompact for BitVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let rv = RawVector::from_compact_serialized(p)?;
        Ok(BitVector::from(rv))
    }
}

impl FromCompact for IntVector {
    fn from_compact_serialized<P: AsRef<Path>>(p: P) -> Result<Self> {
        let mut f = File::open(p)?;
        let _static_flag = u64::read_from(&mut f)?;
        // What is this? this just sets the bits per element
        //assert_eq!(static_flag, 0); // TODO

        let width = u64::read_from(&mut f)? as usize;
        assert!(width > 0);

        let len = u64::read_from(&mut f)? as usize;
        let _capacity = u64::read_from(&mut f)?; // we don't use this

        let words = Vec::read_from(&mut f)?;
        let rv = RawVector::from_parts(len * width, words);

        let iv = IntVector::from_parts(len, width, rv);

        Ok(iv)
    }
}

#[cfg(test)]
mod int_vec_test {
    use super::*;
    use simple_sds::ops::Vector;
    use simple_sds::ops::Access;


    #[test]
    fn test() {
        let pos_fp = "test_data/tiny_index/pos.bin";
        let workdir = env!("CARGO_MANIFEST_DIR");
        let workdir = Path::new(workdir);
        let pos_fp = workdir.join(pos_fp);
        let iv = IntVector::from_compact_serialized(&pos_fp).unwrap();
        let width = get_bits_per_element(&pos_fp);
        assert_eq!(iv.width(), width);
        assert_eq!(iv.width(), 3);
        assert_eq!(iv.len(), 4);

        let correct_ints = vec![3, 2, 0, 1];
        for (loaded, correct) in iv.iter().zip(correct_ints.iter()) {
            assert_eq!(loaded, *correct);
        }
    }
}
#[cfg(test)]
mod read_width_tests {
    use super::*;

    #[test]
    fn test_get_bits_per_element_yeast_chr01() {
        let pos_fp = "test_data/yeast_chr01_index/pos.bin";
        // 221918 kmers in yeastchr01, 221918.log2.ceil()
        test_get_bits_per_element(pos_fp, 18);
    }

    #[test]
    fn test_get_bits_per_element_tiny() {
        let pos_fp = "test_data/tiny_index/pos.bin";
        // Note in line 614 of PufferfishIndexer.cpp:
        //     size_t w = std::log2(tlen) + 1;
        // so 4 kmers in tiny example means, width = 2 + 1
        test_get_bits_per_element(pos_fp, 3);
    }

    fn test_get_bits_per_element(pos_fp: &str, correct_w: usize) {
        let workdir = env!("CARGO_MANIFEST_DIR");
        let workdir = Path::new(workdir);
        let pos_fp = workdir.join(pos_fp);
        let loaded_w = get_bits_per_element(pos_fp);

        assert_eq!(loaded_w, correct_w);
    }
}
