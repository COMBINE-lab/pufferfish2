use serde::de::Deserializer;
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use simple_sds::bit_vector::BitVector;
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Rank, Select, SelectZero, Vector};
use simple_sds::raw_vector::RawVector;
use std::borrow::Borrow;

// Extension traits to serialize simple-sds datastructures since they do not
// derive serde drivers.
// NB / TODO: we could remove this and just use forked simple-sds repo and derive Serde traits there.

pub trait AsSerialize<'a>
where
    Self: Borrow<Self>,
{
    type S: Serialize;
    fn as_serialize(&'a self) -> Self::S;
}

pub trait FromDeserialize<'de> {
    type D: Deserialize<'de>;
    fn from_deserialized(value: Self::D) -> Self;
}

pub fn serialize<'a, S, T>(item: &'a T, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: AsSerialize<'a>,
    T::S: Serialize,
{
    let rvv = item.as_serialize();
    rvv.serialize(serializer)
}

pub fn deserialize<'de, D, T>(deserializer: D) -> Result<T, D::Error>
where
    D: Deserializer<'de>,
    T: FromDeserialize<'de>,
    T::D: Deserialize<'de>,
{
    let rvv = T::D::deserialize(deserializer)?;
    Ok(T::from_deserialized(rvv))
}

// Serde Serializable view structs for simple-sds types
#[derive(Deserialize)]
pub struct RawVecDe {
    pub len: usize,
    pub data: Vec<u64>,
}

#[derive(Serialize)]
pub struct RawVecSer<'rv> {
    pub len: usize,
    pub data: &'rv [u64],
}

impl<'a, T> AsSerialize<'a> for &'a T
where
    T: AsSerialize<'a>,
{
    type S = T::S;
    fn as_serialize(&'a self) -> Self::S {
        (*self).as_serialize()
    }
}

impl<'a> AsSerialize<'a> for RawVector {
    type S = RawVecSer<'a>;
    fn as_serialize(&'a self) -> RawVecSer {
        RawVecSer {
            len: self.len(),
            data: self.as_ref(),
        }
    }
}

#[derive(Serialize)]
pub struct BitVecSer<'a> {
    pub supports_rank: bool,
    pub supports_select: bool,
    pub supports_select_zero: bool,

    #[serde(with = "self")]
    pub data: &'a RawVector,
}

#[derive(Deserialize)]
pub struct BitVecDe {
    pub supports_rank: bool,
    pub supports_select: bool,
    pub supposrts_select_zero: bool,

    #[serde(with = "self")]
    pub data: RawVector,
}
impl<'de> FromDeserialize<'de> for RawVector {
    type D = RawVecDe;
    fn from_deserialized(rvv: Self::D) -> Self {
        Self::from_parts(rvv.len, rvv.data)
    }
}

impl<'a> AsSerialize<'a> for BitVector {
    type S = BitVecSer<'a>;
    fn as_serialize(&'a self) -> BitVecSer<'a> {
        BitVecSer {
            supports_rank: self.supports_rank(),
            supports_select: self.supports_select(),
            supports_select_zero: self.supports_select_zero(),
            data: self.as_ref(),
        }
    }
}

impl FromDeserialize<'_> for BitVector {
    type D = BitVecDe;
    fn from_deserialized(bvv: BitVecDe) -> Self {
        let mut bv = Self::from(bvv.data);
        if bvv.supports_rank {
            bv.enable_rank()
        }
        if bvv.supports_select {
            bv.enable_select()
        }

        if bvv.supposrts_select_zero {
            bv.enable_select_zero()
        }
        bv
    }
}

impl<'a, T> AsSerialize<'a> for Vec<T>
where
    T: 'a,
    T: AsSerialize<'a>,
{
    type S = Vec<<T as AsSerialize<'a>>::S>;
    fn as_serialize(&'a self) -> Self::S {
        self.iter().map(|bv| bv.as_serialize()).collect()
    }
}

impl<'de, T> FromDeserialize<'de> for Vec<T>
where
    T: FromDeserialize<'de>,
{
    type D = Vec<T::D>;
    fn from_deserialized(de: Self::D) -> Self {
        Vec::from_iter(de.into_iter().map(|d| T::from_deserialized(d)))
    }
}

#[derive(Serialize)]
pub struct IntVecSer<'a> {
    pub len: usize,
    pub width: usize,

    #[serde(with = "self")]
    pub data: &'a RawVector,
}

#[derive(Deserialize)]
pub struct IntVecDe {
    pub len: usize,
    pub width: usize,

    #[serde(with = "self")]
    pub data: RawVector,
}

impl<'a> AsSerialize<'a> for IntVector {
    type S = IntVecSer<'a>;
    fn as_serialize(&'a self) -> IntVecSer<'a> {
        IntVecSer {
            len: self.len(),
            width: self.width(),
            data: self.as_ref(),
        }
    }
}

impl FromDeserialize<'_> for IntVector {
    type D = IntVecDe;
    fn from_deserialized(ivv: Self::D) -> Self {
        Self::from_parts(ivv.len, ivv.width, ivv.data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bincode;
    use simple_sds::ops::{Pack, Rank, Resize, Select, SelectZero};

    #[quickcheck]
    fn bv(xs: Vec<u64>, size: usize) -> bool {
        let mut rv = RawVector::from_parts(xs.len() * 64, xs);
        let size = size % 10000;
        rv.resize(size, false);
        let mut bv = BitVector::from(rv);
        bv.enable_rank();
        bv.enable_select();
        bv.enable_select_zero();
        let ser = bincode::serialize(&bv.as_serialize()).unwrap();
        let de = bincode::deserialize(&ser).unwrap();
        let de = BitVector::from_deserialized(de);
        de == bv
    }

    #[quickcheck]
    fn iv(xs: Vec<u64>, size: usize) -> bool {
        let mut iv = IntVector::from(xs);
        let size = size % 10000;
        iv.resize(size, 0);
        iv.pack();
        let ser = bincode::serialize(&iv.as_serialize()).unwrap();
        let de = bincode::deserialize(&ser).unwrap();
        let de = IntVector::from_deserialized(de);
        de == iv
    }
}
