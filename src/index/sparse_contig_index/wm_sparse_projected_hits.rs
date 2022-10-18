// use crate::canonical_kmer::CanonicalKmer;
use super::occs::{EncodedOccs, SampledContigOcc};
use super::wm_ref_walker::{RefWalker, WalkCache};
use super::StreamingCache;
use crate::index::contig::ContigOcc;
use crate::index::{DecodeHit, MappedOrientation, MappedRefPos, PuffQuery, PufferfishBase};
use kmers::naive_impl::Kmer;

#[derive(Debug)]
pub struct SparseProjHits<'a, T> {
    pub index: &'a T,
    // queried kmer
    pub kmer: Kmer,
    // Orientation that the kmer maps to seq
    pub mapped_orientation: MappedOrientation,
    // offset in kmer in canonical orientation of contig
    pub offset: usize,
    pub contig_id: usize,
    pub contig_len: usize,
    pub hits: EncodedOccs<'a>,
}

pub struct DumbIterT<T> {
    _phantom: std::marker::PhantomData<T>,
}
impl<T> Iterator for DumbIterT<T> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        None
    }
}

// https://github.com/rust-lang/rust/issues/41481#issue-223637526
// Workaround derive for Clone requiring all generics to impl Clone
impl<'a, T> Clone for SparseProjHits<'a, T> {
    fn clone(&self) -> Self {
        Self {
            index: self.index,
            kmer: self.kmer.clone(),
            mapped_orientation: self.mapped_orientation,
            offset: self.offset,
            contig_id: self.contig_id,
            contig_len: self.contig_len,
            hits: self.hits.clone(),
        }
    }
}

impl<'a, T> DecodeHit for SparseProjHits<'a, T>
where
    T: PufferfishBase + PuffQuery<'a, HitsT = SparseProjHits<'a, T>>, //     T: Pufferfish,
{
    type IterT = DumbIterT<MappedRefPos>;

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
        assert!(i < self.hits.len());

        let occ = match &self.hits {
            EncodedOccs::Sampled(hits) => self.decode_sampled_contig_occ(hits[i]),
            EncodedOccs::NotSampled(_) => self.decode_nonsamp_contig_occ(i),
            EncodedOccs::SampledWM(hits) => self.decode_sampled_contig_occ(hits.encoded_occs[i]),
            EncodedOccs::NotSampledWM(_) => self.decode_nonsamp_contig_occ(i),
            EncodedOccs::NotEncoded(occs) => occs[i].clone(),
        };

        self.decode_contig_info(&occ)
    }

    fn iter(&self) -> Self::IterT {
        todo!()
    }

    fn as_vec(&self) -> Vec<MappedRefPos> {
        let occs = self.as_contig_occs();
        occs.iter()
            .map(|occ| self.decode_contig_info(occ))
            .collect()
    }
}

#[allow(dead_code)]
impl<'cache, 'a, T> SparseProjHits<'a, T>
where
    T: PufferfishBase + PuffQuery<'a, HitsT = SparseProjHits<'a, T>>,
{
    pub fn as_contig_occs(&self) -> Vec<ContigOcc> {
        match &self.hits {
            EncodedOccs::Sampled(hits) => hits
                .iter()
                .map(|&hit| self.decode_sampled_contig_occ(hit))
                .collect(),
            EncodedOccs::NotSampled(_) => (0..self.len())
                .into_iter()
                .map(|i| self.decode_nonsamp_contig_occ(i))
                .collect(),
            EncodedOccs::SampledWM(hits) => {
                let hits = hits.encoded_occs;
                hits.iter()
                    .map(|&hit| self.decode_sampled_contig_occ(hit))
                    .collect()
            }
            EncodedOccs::NotSampledWM(_) => (0..self.len())
                .into_iter()
                .map(|i| self.decode_nonsamp_contig_occ(i))
                .collect(), // _ => unreachable!()
            EncodedOccs::NotEncoded(occs) => occs.clone(),
        }
    }

    pub fn as_vec_w_cache_cache(
        &self,
        wc: &'cache mut WalkCache<'a, T>,
        sc: &mut StreamingCache,
    ) -> Vec<MappedRefPos> {
        if sc.contig_id == self.contig_id {
            if let Some(occs) = &sc.contig_occs {
                occs.iter()
                    .map(|occ| self.decode_contig_info(occ))
                    .collect()
            } else {
                let occs = self.as_contig_occs_w_cache(wc);
                sc.contig_occs = Some(occs.clone()); // update the cached occs
                occs.iter()
                    .map(|occ| self.decode_contig_info(occ))
                    .collect()
            }
        } else {
            log::warn!("Should never get here?");
            let occs = self.as_contig_occs_w_cache(wc);
            occs.iter()
                .map(|occ| self.decode_contig_info(occ))
                .collect()
        }
    }

    pub fn as_vec_w_cache(&self, cache: &'cache mut WalkCache<'a, T>) -> Vec<MappedRefPos> {
        let occs = self.as_contig_occs_w_cache(cache);
        occs.iter()
            .map(|occ| self.decode_contig_info(occ))
            .collect()
    }

    pub fn as_contig_occs_w_cache(&self, cache: &'cache mut WalkCache<'a, T>) -> Vec<ContigOcc> {
        match &self.hits {
            EncodedOccs::Sampled(hits) => hits
                .iter()
                .map(|&hit| self.decode_sampled_contig_occ(hit))
                .collect(),
            EncodedOccs::NotSampled(_) => (0..self.len())
                .into_iter()
                .map(|i| self.decode_nonsamp_contig_occ_w_cache(i, cache))
                .collect(),
            EncodedOccs::SampledWM(hits) => {
                // No use for pred/succ info, just unwrap
                let hits = hits.encoded_occs;
                hits.iter()
                    .map(|&hit| self.decode_sampled_contig_occ(hit))
                    .collect()
            }
            EncodedOccs::NotSampledWM(_) => (0..self.len())
                .into_iter()
                .map(|i| self.decode_nonsamp_contig_occ_w_cache(i, cache))
                .collect(),
            EncodedOccs::NotEncoded(occs) => occs.clone(),
        }
    }

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

    pub fn decode_nonsamp_contig_occ_w_cache(
        &self,
        occ_rank: usize,
        cache: &'cache mut WalkCache<'a, T>,
    ) -> ContigOcc {
        let mut walker = RefWalker::new_w_cache(self, occ_rank, cache);
        walker.walk_with_cache(cache)
    }

    pub fn decode_sampled_contig_occ(&self, encoded: u64) -> ContigOcc {
        // returns position of the i-th match w.r.t the forward reference
        let occ = SampledContigOcc::from(encoded);
        ContigOcc::from(occ)
    }

    pub fn decode_nonsamp_contig_occ(&self, occ_rank: usize) -> ContigOcc {
        let mut walker = RefWalker::new(self, occ_rank);
        walker.walk()
    }
}
