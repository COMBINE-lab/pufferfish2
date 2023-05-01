// crate::io
use bincode;
use std::io::{Error, ErrorKind, Read, Result};

pub trait SerializeTo<P> {
    // TODO: is bincode::result right?
    fn serialize_to(&self, p: P) -> bincode::Result<()>;
}

pub trait DeserializeFrom<P>
where
    Self: Sized,
{
    fn deserialize_from(p: P) -> bincode::Result<Self>;
}

/******************************************************************************/
// ReadFrom --
//   convenience trait to read bytes into types directly from file
//   pointers
/******************************************************************************/
pub trait ReadFrom {
    fn read_from(f: &mut dyn Read) -> Result<Self>
    where
        Self: Sized;
}

impl ReadFrom for u64 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 8] = [0; 8];
        f.read_exact(&mut val)?;
        let val = u64::from_ne_bytes(val);
        Ok(val)
    }
}

impl ReadFrom for f64 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 8] = [0; 8];
        f.read_exact(&mut val)?;
        let val = f64::from_ne_bytes(val);
        Ok(val)
    }
}

impl ReadFrom for u32 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 4] = [0; 4];
        f.read_exact(&mut val)?;
        let val = u32::from_ne_bytes(val);
        Ok(val)
    }
}

impl ReadFrom for i32 {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        let mut val: [u8; 4] = [0; 4];
        f.read_exact(&mut val)?;
        let val = i32::from_ne_bytes(val);
        Ok(val)
    }
}

impl ReadFrom for Vec<u64> {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        type T = u64;
        let mut buf = Vec::new();
        let n_bytes = f.read_to_end(&mut buf)?;
        let sizeof_t = std::mem::size_of::<T>();

        if (n_bytes % sizeof_t) != 0 {
            let msg = format!(
                "Cannot read {} bytes not divisible by {}",
                n_bytes, sizeof_t
            );
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }
        let n_words = n_bytes / sizeof_t;

        let mut v = Vec::with_capacity(n_words);

        for chunk in buf.chunks_exact(sizeof_t) {
            // unwrap is fine here. 8 bytes guaranteed by chunks_exact
            let chunk = chunk.try_into().unwrap();
            let word = T::from_ne_bytes(chunk);
            v.push(word);
        }
        Ok(v)
    }
}

impl ReadFrom for Vec<u32> {
    fn read_from(f: &mut dyn Read) -> Result<Self> {
        type T = u32;
        let mut buf = Vec::new();
        let n_bytes = f.read_to_end(&mut buf)?;
        let sizeof_t = std::mem::size_of::<T>();

        if (n_bytes % sizeof_t) != 0 {
            let msg = format!(
                "Cannot read {} bytes not divisible by {}",
                n_bytes, sizeof_t
            );
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }
        let n_words = n_bytes / sizeof_t;

        let mut v = Vec::with_capacity(n_words);

        for chunk in buf.chunks_exact(sizeof_t) {
            // unwrap is fine here. 8 bytes guaranteed by chunks_exact
            let chunk = chunk.try_into().unwrap();
            let word = T::from_ne_bytes(chunk);
            v.push(word);
        }
        Ok(v)
    }
}

pub fn read_u64_vec_with_len(f: &mut dyn Read, len: usize) -> Result<Vec<u64>> {
    type T = u64;
    let n_bytes = len * std::mem::size_of::<T>();
    let mut f = f.take(n_bytes as u64);
    Vec::read_from(&mut f)
}

pub fn read_u32_vec_with_len(f: &mut dyn Read, len: usize) -> Result<Vec<u32>> {
    type T = u32;
    let n_bytes = len * std::mem::size_of::<T>();
    let mut f = f.take(n_bytes as u64);
    Vec::read_from(&mut f)
}
