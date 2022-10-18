use bincode;
use delegate::delegate;
use simple_sds::bit_vector::BitVector;
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Access, BitVec, Pack, Rank};
use simple_sds::raw_vector::{AccessRaw, RawVector};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

use super::super::contig::ContigOcc;
use super::super::Pufferfish;
use crate::io::{DeserializeFrom, SerializeTo};
use crate::serde_ext::{AsSerialize, FromDeserialize};

use log::debug;

// mod contig_sampler;
use super::contig_samplers::ContigSamplingStrategy;

use super::occs::{
    EncodedOccs, NotSampledContigOcc, NotSampledWMOcc, SampledContigOcc, SampledWMOcc,
};

use crate::wm::WaveletMatrix;

#[derive(Debug, Eq, PartialEq)]
#[allow(dead_code)]
pub struct SparseContigTable {
    sampled_ctable: Vec<u64>,
    sampled_offsets: IntVector,

    sampled_is_pop_bv: BitVector,
    sampled_succ_wm: WaveletMatrix,
    sampled_wm_offsets: IntVector,

    nonsamp_ctable: Vec<u8>,
    nonsamp_offsets: IntVector,

    nonsamp_is_pop_bv: BitVector,
    nonsamp_pred_wm: WaveletMatrix,
    nonsamp_succ_wm: WaveletMatrix,
    nonsamp_wm_offsets: IntVector,

    is_sampled_bv: BitVector,

    // pop_ctg_thresh: usize,
    ref_names: Vec<String>,
    // _ref_exts: Vec<u32>,
}

pub struct SparseContigTableBuilder<'a, T> {
    index: &'a T,
    sampled_ctable: Vec<u64>,
    sampled_offsets: Vec<usize>,
    sampled_wm_succs: Vec<u8>,
    sampled_wm_offsets: Vec<usize>,

    nonsamp_ctable: Vec<u8>,
    nonsamp_offsets: Vec<usize>,
    nonsamp_wm_succs: Vec<u8>, // (succ, o) pairs
    nonsamp_wm_preds: Vec<u8>, // (pred, o) pairs
    nonsamp_wm_offsets: Vec<usize>,

    pop_sampled_thresh: usize, //threshhold for building O(1) ranksel for sampled contigs
    pop_nonsamp_thresh: usize, //threshhold for building O(1) ranksel for sampled contigs
    sampled_is_pop_bv: RawVector,
    nonsamp_is_pop_bv: RawVector,

    is_sampled_bv: BitVector,
}

#[allow(dead_code)]
impl<'a, T> SparseContigTableBuilder<'a, T>
where
    T: Pufferfish,
{
    delegate! {
        to self.index {
            fn num_contigs(&self) -> usize;
            fn contig_len(&self, contig_id: usize) -> usize;
        }
    }

    pub fn new(
        pi: &'a T,
        strategy: ContigSamplingStrategy,
        pop_sampled_thresh: usize,
        pop_nonsamp_thresh: usize,
    ) -> Self {
        debug!("Getting is_sampled bv");
        let is_sampled_bv = strategy.get_is_sampled_bv(pi);
        let sampled_ctable = Vec::new();
        let sampled_offsets = Vec::new();
        let nonsamp_ctable = Vec::new();
        let nonsamp_offsets = Vec::new();

        let n_sampled = is_sampled_bv.count_ones();
        let n_nonsamp = is_sampled_bv.len() - n_sampled;

        let sampled_is_pop_bv = RawVector::with_len(n_sampled, false);
        let sampled_wm_succs = Vec::new();
        let sampled_wm_offsets = Vec::new();

        let nonsamp_is_pop_bv: RawVector = RawVector::with_len(n_nonsamp, false);
        let nonsamp_wm_succs = Vec::new();
        let nonsamp_wm_preds = Vec::new();
        let nonsamp_wm_offsets = Vec::new();

        Self {
            pop_sampled_thresh,
            pop_nonsamp_thresh,
            index: pi,

            is_sampled_bv,

            sampled_ctable, // encodes succ, o, ref_id, and pos in u64
            sampled_offsets,

            sampled_is_pop_bv,
            sampled_wm_succs,
            sampled_wm_offsets,

            nonsamp_ctable,
            nonsamp_offsets,
            nonsamp_is_pop_bv,
            nonsamp_wm_succs,
            nonsamp_wm_preds,
            nonsamp_wm_offsets,
        }
    }

    #[inline]
    fn get_4bit_succ_o_pair(&self, ctg: &ContigOcc, ctg_len: usize) -> u8 {
        // Get 3-bit encoding of successor
        let hi_bits = self.get_3bit_succ(ctg, ctg_len) as u8;
        let enc = hi_bits << 1;

        if enc == 0b1000 {
            return enc;
        }

        if ctg.is_forward() {
            enc | 1
        } else {
            enc
        }
    }

    #[inline]
    fn get_3bit_pred_o_pair(&self, ctg: &ContigOcc) -> u8 {
        // Get 2-bit encoding of successor
        let hi_bits = self.get_pred(ctg) as u8;
        let enc = hi_bits << 1;
        if ctg.is_forward() {
            enc | 1
        } else {
            enc
        }
    }

    #[inline]
    fn get_3bit_succ(&self, ctg: &ContigOcc, ctg_len: usize) -> u64 {
        // Get 3-bit encoding of successor
        let offset = self.index.get_ref_accum_lens()[ctg.ref_id as usize] as usize;
        let next_ref_start_pos =
            self.index.get_ref_accum_lens()[(ctg.ref_id + 1) as usize] as usize;
        let succ_pos = offset + (ctg.pos as usize) + ctg_len;
        let is_terminal_contig = succ_pos >= next_ref_start_pos;
        if is_terminal_contig {
            0b100
        } else {
            self.get_succ(ctg, ctg_len)
        }
    }

    #[inline]
    fn get_succ(&self, ctg: &ContigOcc, ctg_len: usize) -> u64 {
        // Get 2-bit encoding of successor
        let succ_pos = (ctg.pos as usize) + ctg_len;
        self.index.get_refseq_nuc(ctg.ref_id as usize, succ_pos)
    }

    #[inline]
    fn get_pred(&self, ctg: &ContigOcc) -> u64 {
        // Get 2-bit encoding of predecessor
        let pred_pos = (ctg.pos as usize) - 1;
        self.index.get_refseq_nuc(ctg.ref_id as usize, pred_pos)
    }

    #[inline]
    fn get_curr(&self, ctg: &ContigOcc) -> u64 {
        // Get 2-bit encoding of "curr" nucleotide
        // the nucleotide that matches pred's successor
        // or simply the k-th base on the contig
        // read in the orientaiton of the reference
        let curr_pos = (ctg.pos as usize) + self.index.k() - 1;
        self.index.get_refseq_nuc(ctg.ref_id as usize, curr_pos)
    }

    #[inline]
    fn populate(&mut self) {
        // populates offsets *and* sampled info
        let mut curr_sampled_offset = 0;
        let mut curr_nonsamp_offset = 0;

        let mut curr_sampled_wm_offset = 0;
        let mut curr_nonsamp_wm_offset = 0;
        //TODO save the max of offsets for intvec width, or just let intvec impl pack
        let mut sampled_i = 0;
        let mut nonsamp_i = 0;
        for ctg_id in 0..self.num_contigs() {
            let is_sampled = self.is_sampled_bv.get(ctg_id);
            let slice = self.index.get_encoded_contig_occs(ctg_id).to_vec();
            let decoded_occs: Vec<ContigOcc> = slice.iter().map(|e| ContigOcc::from(*e)).collect();
            let mut encoded_occs = slice.to_vec(); // TODO -> just Vec::new()?
            let num_occs = encoded_occs.len();
            let ctg_len = self.index.contig_len(ctg_id);

            if is_sampled {
                // Update the sampled offset
                self.sampled_offsets.push(curr_sampled_offset);
                curr_sampled_offset += num_occs;

                // Generate and push new encodings
                // TODO this could just be
                // decoded_occs.iter().map(|occ| { ... } ).collect()
                for (i, occ) in decoded_occs.iter().enumerate() {
                    let succ = self.get_3bit_succ(occ, ctg_len);
                    let encoded = SampledContigOcc::from_dense_occ(occ, succ).as_encoded();
                    encoded_occs[i] = encoded
                }
                self.sampled_ctable.append(&mut encoded_occs);

                // Push WM if ctg is popular
                let is_pop = num_occs > self.pop_sampled_thresh;
                if is_pop {
                    self.sampled_wm_offsets.push(curr_sampled_wm_offset);
                    curr_sampled_wm_offset += num_occs;
                    let mut seq: Vec<u8> = decoded_occs
                        .iter()
                        .map(|occ| self.get_4bit_succ_o_pair(occ, ctg_len))
                        .collect();
                    self.sampled_wm_succs.append(&mut seq);
                    self.sampled_is_pop_bv.set_bit(sampled_i, true);
                }
                sampled_i += 1;
            } else {
                // A nonsampled contig is stored either in the WM or as a slice in Vec<u8>
                let is_pop = num_occs > self.pop_nonsamp_thresh;
                // Either put (pred, o) and (succ, o) pair in a WM or in a Vec
                // If in a Vec, then we put a (pred, curr, o) pair.
                if is_pop {
                    // Put in the WM
                    self.nonsamp_wm_offsets.push(curr_nonsamp_wm_offset);
                    curr_nonsamp_wm_offset += num_occs;
                    // The (succ, o) pairs are 4 bits, 3 for succ, one for orientation
                    let mut seq: Vec<u8> = decoded_occs
                        .iter()
                        .map(|occ| self.get_4bit_succ_o_pair(occ, ctg_len))
                        .collect();
                    self.nonsamp_wm_succs.append(&mut seq);

                    // The (pred, o) pairs are 3 bits, 2 for succ, one for orientation
                    let mut seq: Vec<u8> = decoded_occs
                        .iter()
                        .map(|occ| self.get_3bit_pred_o_pair(occ))
                        .collect();
                    self.nonsamp_wm_preds.append(&mut seq);
                    self.nonsamp_is_pop_bv.set_bit(nonsamp_i, true);
                } else {
                    // Put in the Vec
                    self.nonsamp_offsets.push(curr_nonsamp_offset);
                    curr_nonsamp_offset += num_occs;
                    for occ in decoded_occs {
                        let succ = self.get_succ(&occ, ctg_len);
                        let pred = self.get_pred(&occ);
                        let curr = self.get_curr(&occ);
                        let encoded_occ = NotSampledContigOcc {
                            o: occ.orientation,
                            succ,
                            pred,
                            curr,
                        };
                        let encoded = encoded_occ.as_encoded();
                        self.nonsamp_ctable.push(encoded);
                    }
                }
                nonsamp_i += 1;
            }
        }

        // Finish, push the last accumulated offset to all offset vecs
        self.sampled_offsets.push(curr_sampled_offset);
        self.sampled_wm_offsets.push(curr_sampled_wm_offset);
        self.nonsamp_offsets.push(curr_nonsamp_offset);
        self.nonsamp_wm_offsets.push(curr_nonsamp_wm_offset);
    }

    pub fn build(mut self) -> SparseContigTable {
        debug!("Collecting nuc-o pairs");
        self.populate();

        debug!("Packing offsets for non WM supported lists");
        let mut sampled_offsets = IntVector::from(self.sampled_offsets);
        sampled_offsets.pack();
        let mut nonsamp_offsets = IntVector::from(self.nonsamp_offsets);
        nonsamp_offsets.pack();

        debug!("Building sampled_is_pop rank support");
        let mut sampled_is_pop_bv = BitVector::from(self.sampled_is_pop_bv);
        sampled_is_pop_bv.enable_rank();

        // Sampled (succ, o) pairs
        // 8 (succ, o) pairs, plus one sentinel value
        debug!("Building sampled succ WM");
        let sampled_succ_wm = WaveletMatrix::with_ranksel_support(&self.sampled_wm_succs, 9);

        debug!("Packing sampled succ WM offsets");
        let mut sampled_wm_offsets = IntVector::from(self.sampled_wm_offsets);
        sampled_wm_offsets.pack();

        debug!("Building non_samp_is_pop_bv rank support");
        let mut nonsamp_is_pop_bv = BitVector::from(self.nonsamp_is_pop_bv);
        nonsamp_is_pop_bv.enable_rank();

        // Not sampled pairs
        // 8 (succ, o) pairs, plus one sentinel value || 8 (pred, o pairs)
        debug!("Building nonsamp pred WM");
        let nonsamp_pred_wm = WaveletMatrix::with_ranksel_support(&self.nonsamp_wm_preds, 8);
        debug!("Building nonsamp succ WM");
        let nonsamp_succ_wm = WaveletMatrix::with_ranksel_support(&self.nonsamp_wm_succs, 9);

        debug!("Packing nonsamp_wm_offsets");
        let mut nonsamp_wm_offsets = IntVector::from(self.nonsamp_wm_offsets);
        nonsamp_wm_offsets.pack();

        debug!("Returning SparseContigTable");
        SparseContigTable {
            // pop_ctg_thresh: self.pop_ctg_thresh,
            sampled_is_pop_bv,
            sampled_ctable: self.sampled_ctable,
            sampled_offsets,
            sampled_succ_wm,
            sampled_wm_offsets,

            nonsamp_is_pop_bv,
            nonsamp_ctable: self.nonsamp_ctable,
            nonsamp_offsets,
            nonsamp_pred_wm,
            nonsamp_succ_wm,
            nonsamp_wm_offsets,

            is_sampled_bv: self.is_sampled_bv,

            //this probably means that these copied members should be in base
            ref_names: self.index.get_ref_names().to_vec(),
            // _ref_exts: self.index.ctg_table._ref_exts.clone(),
        }
    }
}

impl SparseContigTable {
    pub fn enable_wm_ranksel(&mut self) {
        self.sampled_succ_wm.enable_rank();
        self.sampled_succ_wm.enable_select();
    }

    pub fn num_contigs(&self) -> usize {
        self.is_sampled_bv.len()
    }

    pub fn num_sampled(&self) -> usize {
        self.is_sampled_bv.count_ones()
    }

    pub fn num_not_sampled(&self) -> usize {
        self.is_sampled_bv.len() - self.num_sampled()
    }

    pub fn num_total_occs(&self) -> usize {
        self.num_sampled_occs() + self.num_nonsamp_occs()
    }

    pub fn num_sampled_occs(&self) -> usize {
        self.sampled_ctable.len()
    }

    pub fn num_nonsamp_occs(&self) -> usize {
        let num_wm_nonsamp = self.nonsamp_pred_wm.len();
        let num_not_wm_nonsamp = self.nonsamp_ctable.len();
        num_wm_nonsamp + num_not_wm_nonsamp
    }

    pub fn get_encoded_contig_occs(&self, ctg_id: usize) -> EncodedOccs {
        //TODO
        if self.is_sampled(ctg_id) {
            let ctg_id = self.sampled_id(ctg_id);
            self.get_encoded_sampled_hits(ctg_id)
        } else {
            let ctg_id = self.not_sampled_id(ctg_id);
            self.get_encoded_not_samped_hits(ctg_id)
        }
    }

    // pub fn get_sampled_succs(&self, contig_id: usize) -> &[ {
    //     assert!(self.is_sampled(contig_id));
    //     let id = self.sampled_id(contig_id);
    //     let start = self.sampled_offsets[ctg_id] as usize;
    //     let end = self.sampled_offsets[ctg_id + 1] as usize;
    //     self.sampled_succs[start..end]
    // }

    #[inline]
    pub fn is_sampled(&self, contig_id: usize) -> bool {
        self.is_sampled_bv.get(contig_id)
    }

    // pub fn is_sampled_popular(&self, contig_id: usize) -> bool {
    //     self.sampled_is_pop_bv.get(contig_id)
    // }

    pub fn num_ctg_occs(&self, contig_id: usize) -> usize {
        self.get_encoded_contig_occs(contig_id).len()
    }

    #[inline]
    fn get_encoded_sampled_hits(&self, ctg_id: usize) -> EncodedOccs {
        // let start = self.sampled_offsets[ctg_id] as usize;
        // let end = self.sampled_offsets[ctg_id + 1] as usize;
        let start = self.sampled_offsets.get(ctg_id) as usize;
        let end = self.sampled_offsets.get(ctg_id + 1) as usize;
        if self.sampled_is_pop_bv.get(ctg_id) {
            let ref_pos_words = self.sampled_ctable.get(start..end).unwrap();
            let wm_ctg_id = self.sampled_is_pop_bv.rank(ctg_id);
            let wm_start = self.sampled_wm_offsets.get(wm_ctg_id) as usize;
            let wm_end = self.sampled_wm_offsets.get(wm_ctg_id + 1) as usize;
            let wm_slice = self.sampled_succ_wm.get_slice(wm_start, wm_end);
            let occ = SampledWMOcc {
                encoded_occs: ref_pos_words,
                wm_slice,
            };
            EncodedOccs::SampledWM(occ)
        } else {
            EncodedOccs::Sampled(self.sampled_ctable.get(start..end).unwrap())
        }
    }

    #[inline]
    fn get_encoded_not_samped_hits(&self, ctg_id: usize) -> EncodedOccs {
        if self.nonsamp_is_pop_bv.get(ctg_id) {
            let ctg_id = self.nonsamp_is_pop_bv.rank(ctg_id);
            let start = self.nonsamp_wm_offsets.get(ctg_id) as usize;
            let end = self.nonsamp_wm_offsets.get(ctg_id + 1) as usize;
            let pred_wm_slice = self.nonsamp_pred_wm.get_slice(start, end);
            let succ_wm_slice = self.nonsamp_succ_wm.get_slice(start, end);
            let occ = NotSampledWMOcc {
                pred_wm_slice,
                succ_wm_slice,
            };
            EncodedOccs::NotSampledWM(occ)
        } else {
            let ctg_id = self.nonsamp_is_pop_bv.rank_zero(ctg_id);
            let start = self.nonsamp_offsets.get(ctg_id) as usize;
            let end = self.nonsamp_offsets.get(ctg_id + 1) as usize;
            EncodedOccs::NotSampled(self.nonsamp_ctable.get(start..end).unwrap())
        }
    }

    #[inline]
    fn sampled_id(&self, contig_id: usize) -> usize {
        self.is_sampled_bv.rank(contig_id)
    }

    #[inline]
    fn not_sampled_id(&self, contig_id: usize) -> usize {
        self.is_sampled_bv.rank_zero(contig_id)
    }
}

impl<P: AsRef<Path>> SerializeTo<P> for SparseContigTable {
    fn serialize_to(&self, dir: P) -> bincode::Result<()> {
        let fps = SparseContigTablePaths::new(dir);

        let f = File::create(fps.sampled_occs)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.sampled_ctable)?;

        let f = File::create(fps.sampled_offsets)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.sampled_offsets.as_serialize())?;

        let f = File::create(fps.notsamp_occs)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.nonsamp_ctable)?;

        let f = File::create(fps.notsamp_offsets)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.nonsamp_offsets.as_serialize())?;

        let f = File::create(fps.is_sampled_bv)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.is_sampled_bv.as_serialize())?;

        let f = File::create(fps.ref_names)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.ref_names)?;

        let f = File::create(fps.sampled_is_pop_bv)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.sampled_is_pop_bv.as_serialize())?;

        let f = File::create(fps.sampled_succ_wm)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.sampled_succ_wm)?;

        let f = File::create(fps.sampled_wm_offsets)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.sampled_wm_offsets.as_serialize())?;

        let f = File::create(fps.nonsamp_is_pop_bv)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.nonsamp_is_pop_bv.as_serialize())?;

        let f = File::create(fps.nonsamp_succ_wm)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.nonsamp_succ_wm)?;

        let f = File::create(fps.nonsamp_pred_wm)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.nonsamp_pred_wm)?;

        let f = File::create(fps.nonsamp_wm_offsets)?;
        let w = BufWriter::new(f);
        bincode::serialize_into(w, &self.nonsamp_wm_offsets.as_serialize())?;

        // simple_sds::serialize;
        Ok(())
    }
}

impl<P: AsRef<Path>> DeserializeFrom<P> for SparseContigTable {
    fn deserialize_from(dpath: P) -> bincode::Result<Self> {
        let fps = SparseContigTablePaths::new(dpath);

        let f = File::open(fps.sampled_occs)?;
        let r = BufReader::new(f);
        let sampled_ctable = bincode::deserialize_from(r)?;

        let f = File::open(fps.sampled_offsets)?;
        let r = BufReader::new(f);
        let sampled_offsets = bincode::deserialize_from(r)?;
        let sampled_offsets = IntVector::from_deserialized(sampled_offsets);

        let f = File::open(fps.notsamp_occs)?;
        let r = BufReader::new(f);
        let nonsamp_ctable = bincode::deserialize_from(r)?;

        let f = File::open(fps.notsamp_offsets)?;
        let r = BufReader::new(f);
        let nonsamp_offsets = bincode::deserialize_from(r)?;
        let nonsamp_offsets = IntVector::from_deserialized(nonsamp_offsets);

        let f = File::open(fps.is_sampled_bv)?;
        let r = BufReader::new(f);
        let is_sampled_bv = bincode::deserialize_from(r)?;
        let is_sampled_bv = BitVector::from_deserialized(is_sampled_bv);

        let f = File::open(fps.ref_names)?;
        let r = BufReader::new(f);
        let ref_names = bincode::deserialize_from(r)?;

        let f = File::open(fps.sampled_is_pop_bv)?;
        let r = BufReader::new(f);
        let sampled_is_pop_bv = bincode::deserialize_from(r)?;
        let sampled_is_pop_bv = BitVector::from_deserialized(sampled_is_pop_bv);

        let f = File::open(fps.sampled_succ_wm)?;
        let r = BufReader::new(f);
        let sampled_succ_wm = bincode::deserialize_from(r)?;

        let f = File::open(fps.sampled_wm_offsets)?;
        let r = BufReader::new(f);
        let sampled_wm_offsets = bincode::deserialize_from(r)?;
        let sampled_wm_offsets = IntVector::from_deserialized(sampled_wm_offsets);

        let f = File::open(fps.nonsamp_is_pop_bv)?;
        let r = BufReader::new(f);
        let nonsamp_is_pop_bv = bincode::deserialize_from(r)?;
        let nonsamp_is_pop_bv = BitVector::from_deserialized(nonsamp_is_pop_bv);

        let f = File::open(fps.nonsamp_succ_wm)?;
        let r = BufReader::new(f);
        let nonsamp_succ_wm = bincode::deserialize_from(r)?;

        let f = File::open(fps.nonsamp_pred_wm)?;
        let r = BufReader::new(f);
        let nonsamp_pred_wm = bincode::deserialize_from(r)?;

        let f = File::open(fps.nonsamp_wm_offsets)?;
        let r = BufReader::new(f);
        let nonsamp_wm_offsets = bincode::deserialize_from(r)?;
        let nonsamp_wm_offsets = IntVector::from_deserialized(nonsamp_wm_offsets);

        let index = Self {
            sampled_ctable,
            sampled_offsets,
            sampled_is_pop_bv,
            sampled_succ_wm,
            sampled_wm_offsets,

            nonsamp_ctable,
            nonsamp_offsets,
            nonsamp_is_pop_bv,
            nonsamp_pred_wm,
            nonsamp_succ_wm,
            nonsamp_wm_offsets,
            is_sampled_bv,
            ref_names,
        };
        // simple_sds::serialize;
        Ok(index)
    }
}

struct SparseContigTablePaths {
    _dir: PathBuf,
    sampled_occs: PathBuf,
    sampled_offsets: PathBuf,
    sampled_is_pop_bv: PathBuf,
    sampled_succ_wm: PathBuf,
    sampled_wm_offsets: PathBuf,
    notsamp_occs: PathBuf,
    notsamp_offsets: PathBuf,
    nonsamp_is_pop_bv: PathBuf,
    nonsamp_succ_wm: PathBuf,
    nonsamp_pred_wm: PathBuf,
    nonsamp_wm_offsets: PathBuf,
    is_sampled_bv: PathBuf,
    ref_names: PathBuf,
}

impl SparseContigTablePaths {
    pub fn new<P: AsRef<Path>>(dname: P) -> Self {
        let dir = dname.as_ref();
        Self {
            _dir: PathBuf::from(dir),
            sampled_occs: dir.join(fp::SAMPLED_OCCS),
            sampled_offsets: dir.join(fp::SAMPLED_OFFSETS),

            sampled_is_pop_bv: dir.join(fp::SAMPLED_IS_POP_BV),
            sampled_succ_wm: dir.join(fp::SAMPLED_SUCC_WM),
            sampled_wm_offsets: dir.join(fp::SAMPLED_WM_OFFSETS),

            notsamp_occs: dir.join(fp::NOTSAMP_OCCS),
            notsamp_offsets: dir.join(fp::NOTSAMP_OFFSETS),

            nonsamp_is_pop_bv: dir.join(fp::NOTSAMP_IS_POP_BV),
            nonsamp_pred_wm: dir.join(fp::NOTSAMP_PRED_WM),
            nonsamp_succ_wm: dir.join(fp::NOTSAMP_SUCC_WM),
            nonsamp_wm_offsets: dir.join(fp::NOTSAMP_WM_OFFSETS),

            is_sampled_bv: dir.join(fp::ISSAMP_BV),
            ref_names: dir.join(fp::REF_NAMES),
        }
    }
}

mod fp {
    pub const SAMPLED_OCCS: &str = "sparse_ctg_sampled_occs.bin";
    pub const SAMPLED_OFFSETS: &str = "sparse_ctg_sampled_offsets.bin";
    pub const SAMPLED_IS_POP_BV: &str = "sparse_ctg_sampled_is_pop_bv.bin";
    pub const SAMPLED_SUCC_WM: &str = "sparse_ctg_sampled_succ_wm.bin";
    pub const SAMPLED_WM_OFFSETS: &str = "sprase_ctg_sampled_wm_offsets.bin";

    pub const NOTSAMP_OCCS: &str = "sparse_ctg_notsamp_occs.bin";
    pub const NOTSAMP_OFFSETS: &str = "sparse_ctg_notsamp_offsets.bin";
    pub const NOTSAMP_IS_POP_BV: &str = "sparse_ctg_notsamp_is_pop_bv.bin";
    pub const NOTSAMP_PRED_WM: &str = "sparse_ctg_notsamp_pred_wm.bin";
    pub const NOTSAMP_SUCC_WM: &str = "sparse_ctg_notsamp_succ_wm.bin";
    pub const NOTSAMP_WM_OFFSETS: &str = "sprase_ctg_notsamp_wm_offsets.bin";

    pub const ISSAMP_BV: &str = "sparse_ctg_is_samp_bv.bin";
    pub const REF_NAMES: &str = "ref_names.bin";
}
