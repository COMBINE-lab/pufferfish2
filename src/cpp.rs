// All things to do with Cpp impl of pufferfish
use crate::io::{read_u32_vec_with_len, read_u64_vec_with_len, ReadFrom};
use std::fs::File;
use std::io::{Read, Result};
use std::path::Path;

pub trait FromCereal {
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self>
    where
        Self: Sized;

    fn load_from_cereal_archive<P: AsRef<Path>>(p: P) -> Result<Self>
    where
        Self: Sized,
    {
        let mut f = File::open(p)?;
        Self::read_from_cereal_archive(&mut f)
    }
}

impl FromCereal for Vec<String> {
    // Deserialize Cereal input archive for std::vec of strings
    // <len_of_vec><n_chars><str>...<nchars><str>
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self> {
        let len = u64::read_from(f)? as usize;
        let mut v = Vec::with_capacity(len);
        for _ in 0..len {
            let n_chars = u64::read_from(f)? as usize;
            let mut buf = vec![0u8; n_chars];
            f.read_exact(&mut buf)?;

            let s = String::from_utf8(buf).unwrap();
            v.push(s);
        }
        Ok(v)
    }
}

impl FromCereal for Vec<u32> {
    // Deserialize Cereal input archive for std::vec
    // <len_of_vec><word>...<word>
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self> {
        let len = u64::read_from(f)? as usize;
        read_u32_vec_with_len(f, len)
    }
}

impl FromCereal for Vec<u64> {
    // Deserialize Cereal input archive for std::vec
    // <len_of_vec><word>...<word>
    fn read_from_cereal_archive(f: &mut dyn Read) -> Result<Self> {
        let len = u64::read_from(f)? as usize;
        read_u64_vec_with_len(f, len)
    }
}

pub trait DeserializeFromCpp
where
    Self: Sized,
{
    fn deserialize_from_cpp<P: AsRef<Path>>(p: P) -> Result<Self>;
}
