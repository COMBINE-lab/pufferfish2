use kmers::naive_impl::{Base, Kmer};
use log::debug;
use serde::{Deserialize, Serialize};
use simple_sds::bit_vector::BitVector;
use simple_sds::ops::{BitVec, Rank, Select};
use std::io::Result;
use std::path::Path;

use super::compact::FromCompact;
use super::contig::ContigOrientation;
use super::info::Info;
use super::PufferfishFilePaths;
use crate::boophf::BooPHF;
use crate::cpp::DeserializeFromCpp;
use crate::cpp::FromCereal;
use crate::seq::{SeqVector, SeqVectorSlice};
use crate::serde_ext;

#[derive(Clone, Debug, Serialize, Deserialize, Eq, PartialEq)]
pub struct BaseIndex {
    // Info fields
    // pub info: Info,
    pub index_version: u64,
    pub reference_gfa: Vec<String>,
    // pub sampling_type: PufferfishType,
    pub k: u64,
    pub num_kmers: usize,
    pub num_contigs: usize,
    pub seq_len: usize,
    // have_ref_seq: bool,
    pub have_edge_vec: bool,
    pub seq_hash: String,
    pub name_hash: String,
    pub seq_hash_512: String,
    pub name_hash_512: String,
    pub decoy_seq_hash: String,
    pub decoy_name_hash: String,
    pub num_decoys: usize,
    pub first_decoy_index: usize,
    pub keep_duplicates: bool,

    // Traversal
    pub seq: SeqVector,
    pub last_seq_pos: usize,

    #[serde(skip)]
    pub mphf: BooPHF<u64>,
    #[serde(with = "serde_ext")]
    pub bv: BitVector, //contig boundaries
    pub ref_seq: Option<SeqVector>,
    pub ref_lens: Vec<u32>,
    pub _complete_ref_lens: Vec<u32>,
    pub ref_accum_lens: Vec<u64>,
    // ref_names: Vec<String>,
    // _ref_exts: Vec<u32>,
}

impl BaseIndex {
    pub fn check_query(&self, kmer: &Kmer) {
        if kmer.len() != self.k() {
            panic!(
                "Kmer query for kmer with k={} is not supported. Queries should have k={}",
                kmer.len(),
                self.k()
            )
        }
    }

    pub fn get_refseq(&self, ref_id: usize) -> SeqVectorSlice {
        assert!(self.has_ref_seq());
        assert!(ref_id < self.ref_lens.len());
        let pos = self.ref_accum_lens[ref_id] as usize;
        let len = self.ref_lens[ref_id] as usize;

        self.ref_seq.as_ref().unwrap().slice(pos, pos + len)
    }

    pub fn get_useq_kmer_u64(&self, pos: usize) -> u64 {
        self.seq.get_kmer_u64(pos, self.k as u8)
    }

    pub fn get_useq_kmer(&self, pos: usize) -> Kmer {
        Kmer::from_u64(self.get_useq_kmer_u64(pos), self.k as u8)
    }

    pub fn get_contig_nuc(&self, ctg_id: usize, ctg_o: ContigOrientation, offset: usize) -> Base {
        if let ContigOrientation::Forward = ctg_o {
            let pos = self.contig_start_pos(ctg_id) + offset;
            self.seq.get_base(pos)
        } else {
            // TODO: this fails when I look for the last contig right?
            let pos = self.contig_start_pos(ctg_id + 1) - offset - 1;
            kmers::naive_impl::prelude::complement_base(self.seq.get_base(pos))
        }
    }

    #[allow(dead_code)]
    // #[inline]
    pub fn get_ref_lens(&self) -> &Vec<u32> {
        &self.ref_lens
    }

    pub fn get_ref_accum_lens(&self) -> &[u64] {
        &self.ref_accum_lens[..]
    }

    pub fn has_ref_seq(&self) -> bool {
        self.ref_seq.is_some()
    }

    // fn get_ref_name(&self, id: usize) -> Option<&String> {
    //     self.ref_names.get(id)
    // }

    pub fn k(&self) -> usize {
        self.k as usize
    }

    pub fn num_kmers(&self) -> usize {
        self.num_kmers
    }

    pub fn num_contigs(&self) -> usize {
        self.num_contigs
    }

    pub fn num_refs(&self) -> usize {
        self.ref_lens.len()
    }

    pub fn seq_len(&self) -> usize {
        self.seq_len
    }

    pub fn ref_len(&self, ref_id: usize) -> usize {
        self.ref_lens[ref_id] as usize
    }

    #[inline]
    pub fn contig_start_pos(&self, contig_id: usize) -> usize {
        if contig_id == 0 {
            0
        } else {
            self.bv.select(contig_id - 1).unwrap() as usize + 1
        }
    }

    #[inline]
    pub fn contig_len(&self, contig_id: usize) -> usize {
        // start position of contig on seq
        // let offset = if contig_id == 0 {
        //     0
        // } else {
        //     self.bv.select(contig_id - 1).unwrap() as usize + 1
        // };

        self.contig_start_pos(contig_id + 1) - self.contig_start_pos(contig_id)

        // self.bv.select(contig_id).unwrap() + 1 - offset
    }

    pub fn contig_prefix(&self, contig_id: usize, o: ContigOrientation) -> Kmer {
        // returns the first k-mer of a contig in the given orientation
        let offset = match o {
            ContigOrientation::Forward => 0,
            ContigOrientation::Backward => self.contig_len(contig_id) - (self.k as usize),
        };

        let pos = self.contig_start_pos(contig_id) + offset;

        match o {
            ContigOrientation::Forward => self.get_useq_kmer(pos),
            ContigOrientation::Backward => self.get_useq_kmer(pos).to_reverse_complement(),
        }
    }

    pub fn get_contig_fw_kmer_u64(&self, contig_id: usize, offset: usize) -> u64 {
        // return the "offset"-th kmer of given contig
        // always forward...
        let pos = self.contig_start_pos(contig_id);
        self.get_useq_kmer_u64(pos + offset)
    }
}

impl DeserializeFromCpp for BaseIndex {
    fn deserialize_from_cpp<P: AsRef<Path>>(dir: P) -> Result<Self> {
        let files = PufferfishFilePaths::new(dir);

        let info = Info::load(files.info_json);
        //   let ref_seq = Self::load_ref_seq(&info);
        //   let edge_vec = Self::load_edge_vec(&info);

        debug!("Loading mphf");
        let mphf = BooPHF::<u64>::deserialize_from_cpp(files.mphf)?;

        debug!("Loading seq");
        let seq = SeqVector::from_compact_serialized(files.seq)?;
        assert_eq!(seq.len(), info.seq_len);

        let last_seq_pos = seq.len() - (info.kmer_size as usize);

        debug!("Loading refseq");
        let ref_seq = if info.have_ref_seq {
            Some(SeqVector::from_compact_serialized(files.ref_seq)?)
        } else {
            None
        };

        debug!("Loading reflens");
        //TODO assert that the len of ref_seq is the sum of lens of reference sequences
        let ref_lens = Vec::load_from_cereal_archive(files.ref_lens)?;
        let ref_accum_lens = {
            // TODO: note, we are appending one zero to the front so accessing the offset is easy
            let mut v = Vec::load_from_cereal_archive(files.ref_accum_lens)?;
            v.insert(0, 0);
            v
        };
        let complete_ref_lens = Vec::load_from_cereal_archive(files.complete_ref_lens)?;

        // let (ref_names, ref_exts, ctable) = Self::deserialize_contig_table(files.ctable)?;

        debug!("Loading bv");
        let bv = {
            let mut bv = BitVector::from_compact_serialized(files.rank)?;
            bv.enable_rank();
            bv.enable_select();
            bv
        };

        assert_eq!(bv.len(), info.seq_len);

        debug!("Loaded base index");

        Ok(Self {
            seq,
            last_seq_pos,
            mphf,
            bv,
            ref_seq,
            ref_lens,
            ref_accum_lens, // cumalative length of references
            _complete_ref_lens: complete_ref_lens,

            // Info fields
            index_version: info.index_version,
            reference_gfa: info.reference_gfa,
            k: info.kmer_size,
            num_kmers: info.num_kmers,
            num_contigs: info.num_contigs,
            seq_len: info.seq_len,
            // [omitted] have_ref_seq: info.have_ref_seq,
            // [omitted] num_kmers: usize,
            // [omitted] num_contigs: info.num_contigs,
            // [omitted] seq_len: usize,
            // [omitted] have_ref_seq: info.have_ref_seq,
            have_edge_vec: info.have_edge_vec,
            seq_hash: info.seq_hash,
            name_hash: info.name_hash,
            seq_hash_512: info.seq_hash_512,
            name_hash_512: info.name_hash_512,
            decoy_seq_hash: info.decoy_seq_hash,
            decoy_name_hash: info.decoy_name_hash,
            num_decoys: info.num_decoys,
            first_decoy_index: info.first_decoy_index,
            keep_duplicates: info.keep_duplicates,
        })
    }
}
