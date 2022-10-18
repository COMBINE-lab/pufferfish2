use super::filter_rank_select::{FilterRank, FilterSelect};
use super::occs;
use super::occs::{EncodedOccs, NotSampledContigOcc, SampledContigOcc};
use super::wm_sparse_projected_hits::SparseProjHits;

use crate::index::contig::{ContigOcc, ContigOrientation};
use crate::index::{MappedOrientation, PuffQuery, PufferfishBase};
use kmers::naive_impl::{Base, CanonicalKmer, Kmer};

// use log::trace;

use lru::LruCache;

const NOTSAMP_HANDSHAKE_CURR_MASK: u8 = occs::NOTSAMP_O_MASK | occs::NOTSAMP_PRED_MASK;
const NOTSAMP_HANSHAKE_MASK: u8 = occs::NOTSAMP_O_MASK | occs::NOTSAMP_SUCC_MASK;
const SAMPLED_HANSHAKE_MASK: u64 = occs::SAMPLED_O_MASK | occs::SAMPLED_SUCC_MASK;

#[derive(Debug, Clone)]
enum State {
    Continue {
        offset: usize, // distance travelled
        contig_id: usize,
        // contig_len: usize,
        handshake_rank: usize, // rank of (pred, o) pair
        curr: Base,
        pred: Base,
        o: ContigOrientation, // orientation of curr contig
                              // occs: &'a [u8],
    },
    // Walk only returns ref_id and pos, the orientation is in the (pred, o) or (succ, o) pair
    Finished {
        ref_id: u32,
        pos: u32,
    },
}

#[derive(Debug)]
pub struct RefWalker<'a, T> {
    index: &'a T,
    // curr_occ_rank: usize,
    // curr_contig_id: usize,
    state: State,
    finish_o: ContigOrientation,
    // state: RefWalkerState<'a>,
    // encoded_occs: &'a EncodedOccs<'a>,
    // encoded_occs: &'a [u8],
}

impl<'a, T> RefWalker<'a, T>
where
    T: PufferfishBase + PuffQuery<'a, HitsT = SparseProjHits<'a, T>>,
{
    pub fn new(sph: &SparseProjHits<'a, T>, occ_rank: usize) -> Self {
        // let start_state = Self::get_start_state(sph, occ_rank);

        let finish_o;
        let state = match &sph.hits {
            EncodedOccs::NotSampled(hits) => {
                let occ = NotSampledContigOcc::from(hits[occ_rank]);
                let handshake_rank = hits.filter_rank(occ_rank, NOTSAMP_HANDSHAKE_CURR_MASK);
                finish_o = occ.o;
                // println!(
                //     "ctg_id: {}, o: {:?}, base {}, pred, {}",
                //     sph.contig_id, occ.o, occ.curr, occ.pred
                // );

                State::Continue {
                    offset: 0,
                    // contig_len: sph.contig_len,
                    contig_id: sph.contig_id,
                    handshake_rank,
                    curr: occ.curr,
                    pred: occ.pred,
                    o: occ.o,
                }
            }
            EncodedOccs::NotSampledWM(wm_occs) => {
                let symbol = wm_occs.pred_wm_slice.access(occ_rank);
                let handshake_rank = wm_occs.pred_wm_slice.rank(symbol, occ_rank);
                finish_o = if symbol & 0b1 == 1 {
                    ContigOrientation::Forward
                } else {
                    ContigOrientation::Backward
                };
                let pred = symbol >> 1;

                // Get curr
                let o = finish_o;
                let k = sph.index.k() as usize;
                let curr = sph.index.get_contig_nuc(sph.contig_id, o, k - 1);
                State::Continue {
                    offset: 0,
                    // contig_len: sph.contig_len,
                    contig_id: sph.contig_id,
                    handshake_rank,
                    curr,
                    pred: pred as u64,
                    o,
                }
            }
            _ => panic!("Cannot walk from sampled"),
        };
        RefWalker {
            index: sph.index,
            state,
            finish_o,
        }
    }

    pub fn new_w_cache(
        sph: &SparseProjHits<'a, T>,
        occ_rank: usize,
        cache: &mut WalkCache<'a, T>,
    ) -> Self {
        let finish_o;
        let state = match &sph.hits {
            EncodedOccs::NotSampled(hits) => {
                let occ = NotSampledContigOcc::from(hits[occ_rank]);
                let handshake_rank = hits.filter_rank(occ_rank, NOTSAMP_HANDSHAKE_CURR_MASK);
                finish_o = occ.o;
                // println!(
                //     "ctg_id: {}, o: {:?}, base {}, pred, {}",
                //     sph.contig_id, occ.o, occ.curr, occ.pred
                // );

                State::Continue {
                    offset: 0,
                    // contig_len: sph.contig_len,
                    contig_id: sph.contig_id,
                    handshake_rank,
                    curr: occ.curr,
                    pred: occ.pred,
                    o: occ.o,
                }
            }
            EncodedOccs::NotSampledWM(wm_occs) => {
                let symbol = wm_occs.pred_wm_slice.access(occ_rank);
                let handshake_rank = wm_occs.pred_wm_slice.rank(symbol, occ_rank);
                finish_o = if symbol & 0b1 == 1 {
                    ContigOrientation::Forward
                } else {
                    ContigOrientation::Backward
                };
                let pred = symbol >> 1;

                // Get curr
                let o = finish_o;
                let k = sph.index.k() as usize;
                let curr = cache
                    .get_or_insert_curr(sph.contig_id, o, || {
                        sph.index.get_contig_nuc(sph.contig_id, o, k - 1)
                    })
                    .unwrap();
                State::Continue {
                    offset: 0,
                    // contig_len: sph.contig_len,
                    contig_id: sph.contig_id,
                    handshake_rank,
                    curr,
                    pred: pred as u64,
                    o,
                }
            }
            _ => panic!("Cannot walk from sampled"),
        };
        RefWalker {
            index: sph.index,
            state,
            finish_o,
        }
    }

    #[inline]
    fn is_finished(&self) -> bool {
        // TODO could implement an enum that maintians sampled/not-sampled
        // state with occ *and* encoded occs with the right types
        // self.state.is_finished
        matches!(self.state, State::Finished { .. })
    }

    #[inline]
    fn grp(&self, km: &Kmer) -> SparseProjHits<'a, T> {
        // TODO wtf is this func?
        let km = CanonicalKmer::from(km.clone());
        let hits = self.index.get_ref_pos(&km).unwrap();
        // trace!("grp id: {}", hits.contig_id);
        hits
    }

    fn prev(&mut self) {
        assert!(!self.is_finished());
        let k = self.index.k() as usize;
        if let State::Continue {
            offset, // distance travelled
            contig_id,
            // contig_len,
            handshake_rank, // rank of (pred, o) pair
            curr,
            pred,
            o,
        } = self.state
        {
            let pred_suffix = self.get_handshake_kmer(contig_id, pred, o);
            // prev hits are guaranteed so we don't necessarily need all of get_ref_pos
            let km = CanonicalKmer::from(pred_suffix);
            let pred_hits = self.index.get_ref_pos(&km).unwrap();
            let pred_o = match pred_hits.mapped_orientation {
                MappedOrientation::Forward => ContigOrientation::Forward,
                MappedOrientation::Backward => ContigOrientation::Backward,
            };
            let ctg_len = pred_hits.contig_len;
            let contig_id = pred_hits.contig_id; // NOTE: contig_id, is of the PRED
            match &pred_hits.hits {
                EncodedOccs::Sampled(prev_occs) => {
                    // Get index of the correct occurence of succ-o pair of the predecessor
                    let mask = SAMPLED_HANSHAKE_MASK;
                    let pattern = Self::get_sampled_pred_pattern(curr, pred_o);
                    let occ_rank = prev_occs
                        .filter_select(handshake_rank, mask, pattern)
                        .unwrap();

                    let mut occ = SampledContigOcc::from(prev_occs[occ_rank]);
                    let offset = offset + ctg_len - k + 1;
                    occ.pos += offset as u32;
                    self.state = State::Finished {
                        ref_id: occ.ref_id,
                        pos: occ.pos,
                    };
                }
                EncodedOccs::NotSampled(prev_occs) => {
                    let mask = NOTSAMP_HANSHAKE_MASK;
                    let pattern = Self::get_not_sampled_pred_pattern(curr, pred_o);
                    let occ_rank = prev_occs
                        .filter_select(handshake_rank, mask, pattern)
                        .unwrap();
                    let occ = NotSampledContigOcc::from(prev_occs[occ_rank]);

                    let handshake_rank =
                        prev_occs.filter_rank(occ_rank, NOTSAMP_HANDSHAKE_CURR_MASK);
                    self.state = State::Continue {
                        offset: offset + ctg_len - k + 1,
                        contig_id,
                        handshake_rank,
                        curr: occ.curr,
                        pred: occ.pred,
                        o: occ.o,
                    };
                }
                EncodedOccs::SampledWM(prev_occs) => {
                    let pred_o = if let ContigOrientation::Forward = pred_o {
                        1
                    } else {
                        0
                    };
                    let symbol = ((curr << 1) | pred_o) as u8;
                    let occ_rank = prev_occs.wm_slice.select(symbol, handshake_rank).unwrap();

                    let mut occ = SampledContigOcc::from(prev_occs.encoded_occs[occ_rank]);
                    let offset = offset + ctg_len - k + 1;
                    occ.pos += offset as u32;
                    self.state = State::Finished {
                        ref_id: occ.ref_id,
                        pos: occ.pos,
                    };
                }
                EncodedOccs::NotSampledWM(prev_occs) => {
                    let pred_o = if let ContigOrientation::Forward = pred_o {
                        1
                    } else {
                        0
                    };
                    let symbol = ((curr << 1) | pred_o) as u8;
                    let occ_rank = prev_occs
                        .succ_wm_slice
                        .select(symbol, handshake_rank)
                        .unwrap();

                    // Work on start state of next iteration
                    let symbol = prev_occs.pred_wm_slice.access(occ_rank); // (pred, o) pair of the predecessor
                    let handshake_rank = prev_occs.pred_wm_slice.rank(symbol, occ_rank);
                    let o = if symbol & 0b1 == 1 {
                        ContigOrientation::Forward
                    } else {
                        ContigOrientation::Backward
                    };
                    let pred = symbol >> 1;

                    // Get curr
                    // contig_id: is of the predecessor before the match satement starts
                    let k = self.index.k() as usize;
                    let curr = self.index.get_contig_nuc(contig_id, o, k - 1);
                    self.state = State::Continue {
                        offset: offset + ctg_len - k + 1,
                        contig_id,
                        handshake_rank,
                        curr,
                        pred: pred as u64,
                        o,
                    };
                }
                _ => unimplemented!(),
            }
        }
    }

    fn prev_w_cache(&mut self, cache: &mut WalkCache<'a, T>) {
        assert!(!self.is_finished());
        let k = self.index.k() as usize;
        if let State::Continue {
            offset, // distance travelled
            contig_id,
            // contig_len,
            handshake_rank, // rank of (pred, o) pair
            curr,
            pred,
            o,
        } = self.state
        {
            let pred_suffix = self.get_handshake_kmer_w_cache(contig_id, pred, o, cache);
            // prev hits are guaranteed so we don't necessarily need all of get_ref_pos
            // let km = CanonicalKmer::from(pred_suffix);
            let pred_hits = cache
                .get_or_insert_hit(&pred_suffix, || self.grp(&pred_suffix))
                .unwrap();
            let pred_o = match pred_hits.mapped_orientation {
                MappedOrientation::Forward => ContigOrientation::Forward,
                MappedOrientation::Backward => ContigOrientation::Backward,
            };
            let ctg_len = pred_hits.contig_len;
            let contig_id = pred_hits.contig_id;
            match &pred_hits.hits {
                EncodedOccs::Sampled(prev_occs) => {
                    // Get index of the correct occurence of succ-o pair of the predecessor
                    let mask = SAMPLED_HANSHAKE_MASK;
                    let pattern = Self::get_sampled_pred_pattern(curr, pred_o);
                    let occ_rank = prev_occs
                        .filter_select(handshake_rank, mask, pattern)
                        .unwrap();

                    let mut occ = SampledContigOcc::from(prev_occs[occ_rank]);
                    let offset = offset + ctg_len - k + 1;
                    occ.pos += offset as u32;
                    self.state = State::Finished {
                        ref_id: occ.ref_id,
                        pos: occ.pos,
                    };
                }
                EncodedOccs::NotSampled(prev_occs) => {
                    let mask = NOTSAMP_HANSHAKE_MASK;
                    let pattern = Self::get_not_sampled_pred_pattern(curr, pred_o);
                    let occ_rank = prev_occs
                        .filter_select(handshake_rank, mask, pattern)
                        .unwrap();
                    let occ = NotSampledContigOcc::from(prev_occs[occ_rank]);

                    let handshake_rank =
                        prev_occs.filter_rank(occ_rank, NOTSAMP_HANDSHAKE_CURR_MASK);
                    self.state = State::Continue {
                        offset: offset + ctg_len - k + 1,
                        contig_id,
                        handshake_rank,
                        curr: occ.curr,
                        pred: occ.pred,
                        o: occ.o,
                    };
                }
                EncodedOccs::SampledWM(prev_occs) => {
                    let pred_o = if let ContigOrientation::Forward = pred_o {
                        1
                    } else {
                        0
                    };
                    let symbol = ((curr << 1) | pred_o) as u8;
                    let occ_rank = prev_occs.wm_slice.select(symbol, handshake_rank).unwrap();

                    let mut occ = SampledContigOcc::from(prev_occs.encoded_occs[occ_rank]);
                    let offset = offset + ctg_len - k + 1;
                    occ.pos += offset as u32;
                    self.state = State::Finished {
                        ref_id: occ.ref_id,
                        pos: occ.pos,
                    };
                }
                EncodedOccs::NotSampledWM(prev_occs) => {
                    let pred_o = if let ContigOrientation::Forward = pred_o {
                        1
                    } else {
                        0
                    };
                    let symbol = ((curr << 1) | pred_o) as u8;
                    let occ_rank = prev_occs
                        .succ_wm_slice
                        .select(symbol, handshake_rank)
                        .unwrap();

                    // Work on start state of next iteration
                    let symbol = prev_occs.pred_wm_slice.access(occ_rank); // (pred, o) pair of the predecessor
                    let handshake_rank = prev_occs.pred_wm_slice.rank(symbol, occ_rank);
                    let o = if symbol & 0b1 == 1 {
                        ContigOrientation::Forward
                    } else {
                        ContigOrientation::Backward
                    };
                    let pred = symbol >> 1;

                    // Get curr
                    let k = self.index.k() as usize;
                    let curr = cache
                        .get_or_insert_curr(contig_id, o, || {
                            self.index.get_contig_nuc(contig_id, o, k - 1)
                        })
                        .unwrap();
                    self.state = State::Continue {
                        offset: offset + ctg_len - k + 1,
                        contig_id,
                        handshake_rank,
                        curr,
                        pred: pred as u64,
                        o,
                    };
                }
                _ => unimplemented!(),
            }
        }
    }

    #[inline]
    fn get_sampled_pred_pattern(curr: Base, pred_o: ContigOrientation) -> u64 {
        // TODO: could make this "faster" by just manipulating bits
        // but instantiating the pred occ with params is ok for now
        let pred_occ = SampledContigOcc {
            succ: curr,
            o: pred_o,
            ref_id: 0,
            pos: 0,
        };
        pred_occ.as_encoded()
    }

    #[inline]
    fn get_not_sampled_pred_pattern(curr: Base, pred_o: ContigOrientation) -> u8 {
        // TODO: could make this "faster" by just manipulating bits (maybe not though...)
        // but instantiating the pred occ with params is more correct
        let pred_occ = NotSampledContigOcc {
            succ: curr,
            o: pred_o,
            pred: 0,
            curr: 0,
        };
        pred_occ.as_encoded()
    }

    #[inline]
    fn get_handshake_kmer_w_cache(
        &self,
        ctg_id: usize,
        pred: Base,
        o: ContigOrientation,
        cache: &mut WalkCache<'a, T>,
    ) -> Kmer {
        let km = cache.get_or_insert_prefix(ctg_id, o, pred, || {
            let mut km = self.index.contig_prefix(ctg_id, o);
            km.prepend_base(pred);
            km
        });
        km.unwrap().clone() //TODO why unwrap?
    }

    #[inline]
    fn get_handshake_kmer(&self, ctg_id: usize, pred: Base, o: ContigOrientation) -> Kmer {
        let mut km = self.index.contig_prefix(ctg_id, o);
        km.prepend_base(pred);
        km
    }

    pub fn walk(&mut self) -> ContigOcc {
        loop {
            if let State::Finished { ref_id, pos } = self.state {
                break ContigOcc {
                    ref_id,
                    pos,
                    orientation: self.finish_o,
                };
            } else {
                self.prev();
            }
        }
    }

    pub fn walk_with_cache(&mut self, cache: &mut WalkCache<'a, T>) -> ContigOcc {
        loop {
            if let State::Finished { ref_id, pos } = self.state {
                break ContigOcc {
                    ref_id,
                    pos,
                    orientation: self.finish_o,
                };
            } else {
                self.prev_w_cache(cache);
            }
        }
    }

}

pub struct WalkCache<'a, T> {
    // (contig ID, Orientation, Predecessor)
    prefixes: LruCache<(usize, ContigOrientation, Base), Kmer>,
    curr_bases: LruCache<(usize, ContigOrientation), Base>,
    hits: LruCache<u64, SparseProjHits<'a, T>>,
}

impl<'a, T> WalkCache<'a, T> {
    pub fn new(n: usize) -> Self {
        Self {
            prefixes: LruCache::new(n),
            curr_bases: LruCache::new(n),
            hits: LruCache::new(n),
        }
    }

    pub fn get_or_insert_curr<F>(
        &mut self,
        contig_id: usize,
        o: ContigOrientation,
        f: F,
    ) -> Option<Base>
    where
        F: Fn() -> Base,
    {
        self.curr_bases.get_or_insert((contig_id, o), f).cloned()
    }

    pub fn get_or_insert_prefix<F>(
        &mut self,
        contig_id: usize,
        o: ContigOrientation,
        p: Base,
        f: F,
    ) -> Option<&Kmer>
    where
        F: Fn() -> Kmer,
    {
        self.prefixes.get_or_insert((contig_id, o, p), f)
    }

    pub fn get_or_insert_hit<'s, F>(
        &'s mut self,
        fw_km: &Kmer,
        f: F,
    ) -> Option<&'s SparseProjHits<'a, T>>
    where
        F: Fn() -> SparseProjHits<'a, T>,
    {
        self.hits.get_or_insert(fw_km.into_u64(), f)
    }
}
