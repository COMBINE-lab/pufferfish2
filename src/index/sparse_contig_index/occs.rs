use crate::index::contig::{ContigOcc, ContigOrientation};
use crate::wm::WaveletMatrixSlice;
use kmers::naive_impl::Base;

// Bitwise layout of not NOT SAMPLED u8 encoded contig occurence from low to high order bits
//  orientation := 1 bit,
//  pred        := 2 bits
//  succ        := 2 bits
//  curr        := 2 bits
// for convienience, the "k"-th nucleotide read on the reference is stored as the
// "curr" base so as to avoid lookups in useq
#[allow(dead_code)]
pub const NOTSAMP_O_OFFSET: usize = 0;
pub const NOTSAMP_PRED_OFFSET: usize = 1;
pub const NOTSAMP_SUCC_OFFSET: usize = 3;
pub const NOTSAMP_CURR_OFFSET: usize = 5;

pub const NOTSAMP_O_MASK: u8 = 0b00000001;
pub const NOTSAMP_PRED_MASK: u8 = 0b00000110;
pub const NOTSAMP_SUCC_MASK: u8 = 0b00011000;
pub const NOTSAMP_CURR_MASK: u8 = 0b01100000;
#[allow(dead_code)]
pub const NOTSAMP_MASK: u8 = 0b01111111;
// TODO: We can mask before decoding a notsampled contig occurence.
//       But correct implementations should not care about 7th bit.

// Bitwise layout of SAMPLED u64 encoded contig occurences from low to high order bits
//  orientation := 1 bit,
//  succ        := 3 bits,
//  ref_id      := 28 bits,
//  pos         := 32 bits,
// Not all successors are valid, as some Sampled contigs have no successors.
// The three bit successor is valid iff the higheset order bit is zero.
// 0b100 indicates a "null" successor
#[allow(dead_code)]
pub const SAMPLED_O_OFFSET: usize = 0;
pub const SAMPLED_SUCC_OFFSET: usize = 1;
pub const SAMPLED_REF_ID_OFFSET: usize = 4;
pub const SAMPLED_POS_OFFSET: usize = 32;

#[allow(dead_code)]
pub const SAMPLED_O_MASK: u64 = 0b1;
pub const SAMPLED_SUCC_MASK: u64 = 0b1110;
pub const SAMPLED_REF_ID_MASK: u64 = 0xFFFFFFF0;
pub const SAMPLED_POS_MASK: u64 = 0xFFFFFFFF00000000;

#[derive(Debug, Clone)]
pub enum EncodedOccs<'a> {
    NotEncoded(Vec<ContigOcc>),
    Sampled(&'a [u64]),
    NotSampled(&'a [u8]),
    SampledWM(SampledWMOcc<'a>),
    NotSampledWM(NotSampledWMOcc<'a>),
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct SampledWMOcc<'a> {
    // pub succ_wm: &'a WaveletMatrix, // pointer into wavelet matrix
    // pub wm_s: usize, // start position in wavelet matrix
    // pub wm_e: usize, // end position in the wavelt matrix
    pub wm_slice: WaveletMatrixSlice<'a>,
    pub encoded_occs: &'a [u64], //
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct NotSampledWMOcc<'a> {
    // pub succ_wm: &'a WaveletMatrix, // pointer into wavelet matrix
    // pub wm_s: usize, // start position in wavelet matrix
    // pub wm_e: usize, // end position in the wavelt matrix
    pub pred_wm_slice: WaveletMatrixSlice<'a>,
    pub succ_wm_slice: WaveletMatrixSlice<'a>,
    // pub encoded_occs: &'a [u64], //
}

impl SampledWMOcc<'_> {
    pub fn len(&self) -> usize {
        self.wm_slice.len()
    }

    pub fn is_empty(&self) -> bool {
        self.wm_slice.is_empty()
    }

    pub fn access_succ_o(&self, i: usize) -> u8 {
        self.wm_slice.access(i)
    }

    //     pub fn select_succ_o(&self, a: u8, rank: usize) -> Option<usize> {
    //         self.wm_slice.select(a, rank)
    //     }
}

#[derive(Debug, Eq, PartialEq)]
pub struct SampledContigOcc {
    pub succ: Base,
    pub o: ContigOrientation,
    pub pos: u32,
    pub ref_id: u32,
}

#[derive(Debug, Eq, PartialEq)]
pub struct NotSampledContigOcc {
    pub pred: Base,
    pub succ: Base,
    pub o: ContigOrientation,
    pub curr: Base,
}

impl EncodedOccs<'_> {
    pub fn len(&self) -> usize {
        match self {
            Self::NotEncoded(hits) => hits.len(),
            Self::Sampled(hits) => hits.len(),
            Self::NotSampled(hits) => hits.len(),
            Self::SampledWM(hits_wm) => hits_wm.encoded_occs.len(),
            Self::NotSampledWM(hits_wm) => hits_wm.succ_wm_slice.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    // pub fn is_sampled(&self) -> bool {
    //     matches!(*self, EncodedOccs::Sampled(_))
    // }

    // pub fn is_sampled_wm(&self) -> bool {
    //     matches!(*self, EncodedOccs::SampledWM(_))
    // }

    // pub fn is_not_sampled(&self) -> bool {
    //     matches!(*self, EncodedOccs::NotSampled(_))
    // }

    // pub fn is_not_sampled_wm(&self) -> bool {

    // }

    pub fn unwrap_sampled(&self) -> &[u64] {
        match self {
            EncodedOccs::Sampled(occs) => occs,
            _ => panic!("Cannot unwrap_sampled not Encoded::Sampled"),
        }
    }

    pub fn unwrap_sampled_wm(&self) -> &SampledWMOcc {
        match self {
            EncodedOccs::SampledWM(occs) => occs,
            _ => panic!("Cannot unwrap_sampled not EncodedOccs::SampledWM"),
        }
    }

    pub fn unwrap_not_sampled(&self) -> &[u8] {
        match self {
            EncodedOccs::NotSampled(occs) => occs,
            _ => panic!("Cannot unwrap_not_sampled not Encoded::NotSampled"),
        }
    }
}

impl From<SampledContigOcc> for ContigOcc {
    fn from(occ: SampledContigOcc) -> Self {
        Self {
            ref_id: occ.ref_id as u32,
            pos: occ.pos as u32,
            orientation: occ.o,
        }
    }
}

impl SampledContigOcc {
    pub fn from_dense_occ(occ: &ContigOcc, succ: u64) -> Self {
        Self {
            succ,
            pos: occ.pos,
            ref_id: occ.ref_id,
            o: occ.orientation,
        }
    }

    pub fn as_encoded(&self) -> u64 {
        let o = self.o.to_u64();
        let succ = self.succ << SAMPLED_SUCC_OFFSET;
        let ref_id = (self.ref_id << SAMPLED_REF_ID_OFFSET) as u64;
        let pos = (self.pos as u64) << SAMPLED_POS_OFFSET;
        pos | ref_id | succ | o
    }
}

impl From<u64> for SampledContigOcc {
    fn from(encoded: u64) -> Self {
        // Note, you could shift then mask, but mask-then-shift fails faster if there are bugs
        let succ = (encoded & SAMPLED_SUCC_MASK) >> SAMPLED_SUCC_OFFSET;
        let ref_id = ((encoded & SAMPLED_REF_ID_MASK) >> SAMPLED_REF_ID_OFFSET) as u32;
        let pos = ((encoded & SAMPLED_POS_MASK) >> SAMPLED_POS_OFFSET) as u32;
        let o = ContigOrientation::from_u64(encoded & SAMPLED_O_MASK);
        Self {
            o,
            succ,
            ref_id,
            pos,
        }
    }
}

impl NotSampledContigOcc {
    pub fn as_encoded(&self) -> u8 {
        let pred = (self.pred << NOTSAMP_PRED_OFFSET) as u8;
        let succ = (self.succ << NOTSAMP_SUCC_OFFSET) as u8;
        let curr = (self.curr << NOTSAMP_CURR_OFFSET) as u8;
        let o = self.o.to_u8();
        pred | succ | curr | o
    }
}

impl From<u8> for NotSampledContigOcc {
    fn from(encoded: u8) -> Self {
        // // TODO: masking off unneeded bits?
        // let encoded = encoded & NOTSAMP_MASK;
        // Note, you could shift then mask, but masking then shifting fails faster if there are bugs
        let pred = ((encoded & NOTSAMP_PRED_MASK) >> NOTSAMP_PRED_OFFSET) as Base;
        let succ = ((encoded & NOTSAMP_SUCC_MASK) >> NOTSAMP_SUCC_OFFSET) as Base;
        let curr = ((encoded & NOTSAMP_CURR_MASK) >> NOTSAMP_CURR_OFFSET) as Base;
        let o = ContigOrientation::from_u8(encoded & NOTSAMP_O_MASK);
        Self {
            pred,
            succ,
            curr,
            o,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consts() {
        assert_eq!(NOTSAMP_O_MASK, 0b1 << NOTSAMP_O_OFFSET);
        assert_eq!(NOTSAMP_PRED_MASK, 0b11 << NOTSAMP_PRED_OFFSET);
        assert_eq!(NOTSAMP_SUCC_MASK, 0b11 << NOTSAMP_SUCC_OFFSET);
        assert_eq!(NOTSAMP_CURR_MASK, 0b11 << NOTSAMP_CURR_OFFSET);

        assert_eq!(
            NOTSAMP_O_MASK ^ NOTSAMP_PRED_MASK ^ NOTSAMP_SUCC_MASK ^ NOTSAMP_CURR_MASK,
            NOTSAMP_MASK
        );

        assert_eq!(SAMPLED_O_MASK, 0b1 << SAMPLED_O_OFFSET);
        assert_eq!(SAMPLED_SUCC_MASK, 0b111 << SAMPLED_SUCC_OFFSET);
        assert_eq!(SAMPLED_REF_ID_MASK, 0x0FFFFFFF << SAMPLED_REF_ID_OFFSET);
        assert_eq!(SAMPLED_POS_MASK, 0xFFFFFFFF << SAMPLED_POS_OFFSET);

        assert_eq!(
            SAMPLED_O_MASK ^ SAMPLED_SUCC_MASK ^ SAMPLED_REF_ID_MASK ^ SAMPLED_POS_MASK,
            u64::MAX
        );
    }

    #[quickcheck]
    fn encode_decode_sampled_identity(encoded: u64) -> bool {
        SampledContigOcc::from(encoded).as_encoded() == encoded
    }

    #[quickcheck]
    fn encode_decode_not_sampled_identity(encoded: u8) -> bool {
        let encoded = encoded & NOTSAMP_MASK;
        NotSampledContigOcc::from(encoded).as_encoded() == encoded
    }

    #[test]
    fn not_sampled_encoding() {
        let encoding = 0b00110110;
        let decoded = NotSampledContigOcc::from(encoding);
        let ans = NotSampledContigOcc {
            o: ContigOrientation::Backward,
            pred: 0b11,
            succ: 0b10,
            curr: 0b01,
        };

        assert_eq!(decoded, ans);
        assert_eq!(encoding, ans.as_encoded());
    }

    #[test]
    fn sampled_encoding() {
        // ENCODED 0x <31 63 1c 94> <7a a4 8c c> 5
        // last four bits: 010 1
        let encoding = 0x31631c947aa48cc5;
        let decoded = SampledContigOcc::from(encoding);
        let ans = SampledContigOcc {
            o: ContigOrientation::Forward,
            succ: 0b10, //C
            ref_id: 0x7aa48cc,
            pos: 0x31631c94,
        };

        assert_eq!(decoded, ans);
        assert_eq!(encoding, ans.as_encoded());
    }
}
