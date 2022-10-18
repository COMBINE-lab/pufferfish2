#[derive(Debug, Eq, PartialEq, Clone, Copy, Hash)]
pub enum ContigOrientation {
    Forward,
    Backward,
}

impl ContigOrientation {
    // Could also use num_derive crate or impl FromPrimitive...
    pub fn to_u8(self) -> u8 {
        match self {
            ContigOrientation::Forward => 1,
            ContigOrientation::Backward => 0,
        }
    }

    pub fn to_u64(self) -> u64 {
        match self {
            ContigOrientation::Forward => 1,
            ContigOrientation::Backward => 0,
        }
    }

    pub fn from_u8(w: u8) -> Self {
        if w > 0 {
            ContigOrientation::Forward
        } else {
            ContigOrientation::Backward
        }
    }

    pub fn from_u64(w: u64) -> Self {
        if w > 0 {
            ContigOrientation::Forward
        } else {
            ContigOrientation::Backward
        }
    }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct ContigOcc {
    // Bitwise layout of SAMPLED u64 encoded contig occurences from low to high order bits
    // ref_id := 32 bits
    // pos    := 31 bits
    pub ref_id: u32,
    pub pos: u32,
    pub orientation: ContigOrientation,
}

impl ContigOcc {
    pub fn is_forward(&self) -> bool {
        self.orientation == ContigOrientation::Forward
    }

    pub fn is_backward(&self) -> bool {
        !self.is_forward()
    }
}

impl From<u64> for ContigOcc {
    fn from(word: u64) -> Self {
        let ref_id = (word & 0xFFFFFFFF) as u32;
        let pos_word = word >> 32;
        let pos = (pos_word & 0x7FFFFFFF) as u32;
        let orientation = if (pos_word & 0x80000000) > 0 {
            ContigOrientation::Forward
        } else {
            ContigOrientation::Backward
        };
        Self {
            ref_id,
            pos,
            orientation,
        }
    }
}
