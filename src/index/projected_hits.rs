// #![allow(incomplete_features)]
// #![feature(generic_associated_types)]

use super::contig::ContigOcc;
use super::{DecodeHit, MappedOrientation, MappedRefPos};
use kmers::naive_impl::Kmer;
// use log::warn;

#[derive(Debug, Clone)]
pub struct ProjectedHits<'a> {
    // queried kmer
    pub kmer: Kmer,
    // Orientation that the kmer maps to seq
    pub mapped_orientation: MappedOrientation,
    // offset in kmer in canonical orientation of contig
    pub offset: usize,
    pub contig_id: usize,
    pub contig_len: usize,
    pub hits: &'a [u64],
}

pub struct ProjectedHitsIterator<'a> {
    curr: usize,
    projected_hits: ProjectedHits<'a>,
    // projected_hits: &'a ProjectedHits<'a> // Streaming iterators are hard in rust.
}
impl<'a> DecodeHit for ProjectedHits<'a> {
    type IterT = ProjectedHitsIterator<'a>;

    fn len(&self) -> usize {
        self.hits.len()
    }

    fn contig_len(&self) -> usize {
        self.contig_len
    }

    fn contig_id(&self) -> usize {
        self.contig_id
    }

    fn mapped_orientation(&self) -> MappedOrientation {
        self.mapped_orientation
    }

    fn is_empty(&self) -> bool {
        // TODO: kind of silly as this should never be empty
        self.hits.is_empty()
    }

    fn decode_hit(&self, i: usize) -> MappedRefPos {
        // returns position of the i-th match w.r.t the forward reference
        assert!(i < self.hits.len());
        self.decode_encoded(self.hits[i])
    }
    fn iter(&self) -> Self::IterT {
        todo!()
        // warn!("For faster iteration without early stopping, use as_vec(...) instead!");
        // ProjectedHitsIterator {
        //     projected_hits: self.clone(),
        //     curr: 0,
        // }
    }

    fn as_vec(&self) -> Vec<MappedRefPos> {
        self.hits
            .iter()
            .map(|&hit| self.decode_encoded(hit))
            .collect()
    }
}

impl ProjectedHits<'_> {
    #[inline]
    pub fn decode_encoded(&self, encoded: u64) -> MappedRefPos {
        self.decode_contig_info(&ContigOcc::from(encoded))
    }

    #[inline]
    pub fn decode_contig_info(&self, contig_info: &ContigOcc) -> MappedRefPos {
        let ref_id = contig_info.ref_id as usize;
        let k = self.kmer.k as usize;
        let contig_pos = contig_info.pos as usize;

        let pos = if contig_info.is_forward() {
            self.offset + contig_pos
        } else {
            contig_pos + (self.contig_len - self.offset) - k
        };

        let orientation = if contig_info.is_forward() {
            self.mapped_orientation
        } else {
            self.mapped_orientation.reverse()
        };

        MappedRefPos {
            ref_id,
            pos,
            orientation,
        }
    }
}

// impl<'a, 'b> ProjectedHitsIterator<'a> {
//     fn new(projected_hits: &'a ProjectedHits<'a>) -> Self {
//         Self {
//             projected_hits,
//             curr: 0,
//         }
//     }
// }

// impl<'a, 'b> IntoIterator for &'a ProjectedHits<'a> {
//     type Item = MappedRefPos;
//     type IntoIter = ProjectedHitsIterator<'a>;
//     fn into_iter(self) -> Self::IntoIter {
//         self.iter()
//     }
// }

impl<'a> Iterator for ProjectedHitsIterator<'a> {
    type Item = MappedRefPos;
    fn next(&mut self) -> Option<Self::Item> {
        if self.curr < self.projected_hits.len() {
            let mapped_pos = Some(self.projected_hits.decode_hit(self.curr));
            self.curr += 1;
            mapped_pos
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpp::DeserializeFromCpp;
    use crate::index::{CanonicalKmer, DenseIndex, PuffQuery};
    use crate::test_utils::*;

    #[test]
    fn tiny_decode_one_hit() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        // --ACC-
        // AAACCC
        let mut acc = CanonicalKmer::from("acc");
        let acc_hits = pi.get_ref_pos(&acc).unwrap();
        let rp = MappedRefPos::new(0, 2, MappedOrientation::Forward);

        assert_eq!(acc_hits.len(), 1);
        assert_eq!(acc_hits.decode_hit(0), rp);

        acc.swap();
        let acc_hits = pi.get_ref_pos(&acc).unwrap();
        let rp = MappedRefPos::new(0, 2, MappedOrientation::Backward);

        assert_eq!(acc_hits.len(), 1);
        assert_eq!(acc_hits.decode_hit(0), rp);
    }

    #[test]
    fn tiny_rc_decode_one_hit() {
        let p = to_abs_path(TINY_RC_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        // -CCA--
        // GGGTTT
        let mut acc = CanonicalKmer::from("acc");
        let acc_hits = pi.get_ref_pos(&acc).unwrap();
        let rp = MappedRefPos::new(0, 1, MappedOrientation::Backward);

        assert_eq!(acc_hits.len(), 1);
        assert_eq!(acc_hits.decode_hit(0), rp);

        acc.swap();
        let acc_hits = pi.get_ref_pos(&acc).unwrap();
        let rp = MappedRefPos::new(0, 1, MappedOrientation::Forward);

        assert_eq!(acc_hits.len(), 1);
        assert_eq!(acc_hits.decode_hit(0), rp);
    }

    #[test]
    fn yeast_decode_one_hit() {
        let p = to_abs_path(YEAST_CHR01_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        // --ACC-
        // AAACCC
        let seq = "ccacaccacacccacacacccacacaccaca";
        let kmer = CanonicalKmer::from(seq);
        let hits = pi.get_ref_pos(&kmer).unwrap();
        let rp = MappedRefPos::new(0, 0, MappedOrientation::Forward);
        assert_eq!(hits.decode_hit(0), rp);
    }

    #[test]
    fn yeast_decode_one_hit_from_as_vec() {
        let p = to_abs_path(YEAST_CHR01_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        // --ACC-
        // AAACCC
        let seq = "ccacaccacacccacacacccacacaccaca";
        let kmer = CanonicalKmer::from(seq);
        let hits = pi.get_ref_pos(&kmer).unwrap();
        let rp = MappedRefPos::new(0, 0, MappedOrientation::Forward);
        let mut ctr = 0;
        for ip in hits.as_vec() {
            assert_eq!(ip, rp);
            ctr += 1;
        }
        assert_eq!(ctr, 1);
    }
}

impl<'a, OtherT> PartialEq<OtherT> for ProjectedHits<'a>
where
    OtherT: DecodeHit,
{
    fn eq(&self, other: &OtherT) -> bool {
        self.as_vec() == other.as_vec()
    }
}
