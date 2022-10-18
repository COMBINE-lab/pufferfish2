use delegate::delegate;
use log::debug;
use simple_sds::int_vector::IntVector;
use simple_sds::ops::{Access, Rank, Vector};
use std::io::Result;
use std::path::Path;

use super::base_index::BaseIndex;
use super::compact::FromCompact;

use super::info::Info;
use super::projected_hits::ProjectedHits;
use super::refseq_iter::RefIterator;
use super::{
    CachedQuery, DenseContigTable, MappedOrientation, PuffQuery, PufferfishBase,
    PufferfishFilePaths, PufferfishType, QueryCache,
};
use crate::boophf::MPHF;
use crate::cpp::DeserializeFromCpp;
use kmers::naive_impl::{CanonicalKmer, MatchType};

#[derive(Clone)]
pub struct DenseIndex {
    // all fields are pub(super) for easy unpacking in super/crate::index
    pub(super) pos: IntVector,
    pub(super) base: BaseIndex,
    pub(super) ctg_table: DenseContigTable,
}

impl DenseIndex {
    fn get_kmer_pos(&self, kmer: &CanonicalKmer) -> Option<(usize, MappedOrientation)> {
        // return the position of queried kmer
        // on seq, concatenated sequence of contigs
        self.check_query(&kmer.get_fw_mer());
        let canonical_query = kmer.get_canonical_kmer();

        let word: u64 = canonical_query.into_u64();
        let hash = self.base.mphf.lookup(&word);

        match hash {
            None => None,
            Some(i) => {
                let candidate_pos = self.pos.get(i as usize) as usize;
                let kmer_at_pos = self
                    .base
                    .seq
                    .get_kmer(candidate_pos as usize, self.base.k as u8);

                match kmer.get_kmer_equivalency(&kmer_at_pos) {
                    MatchType::IdentityMatch => Some((candidate_pos, MappedOrientation::Forward)),
                    MatchType::TwinMatch => Some((candidate_pos, MappedOrientation::Backward)),
                    _ => None,
                }
            }
        }
    }
}

impl DeserializeFromCpp for DenseIndex {
    fn deserialize_from_cpp<P: AsRef<Path>>(dir: P) -> Result<Self> {
        debug!("Loading base index");
        if let Ok(base_struct) = BaseIndex::deserialize_from_cpp(&dir) {
            let files = PufferfishFilePaths::new(&dir);

            debug!("Loading pos");
            let pos = IntVector::from_compact_serialized(files.pos)?;
            assert_eq!(pos.len(), base_struct.num_kmers);

            debug!("Loading contig table");
            let ctg_table = DenseContigTable::deserialize_from_cpp(&dir)?;

            debug!("Loaded dense index");

            Ok(Self {
                pos,
                ctg_table,
                base: base_struct,
            })
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "could not load index.",
            ))
        }
    }
}

impl PufferfishBase for DenseIndex {
    fn num_total_occs(&self) -> usize {
        self.ctg_table.num_total_occs()
    }
    fn num_ctg_occs(&self, contig_id: usize) -> usize {
        self.ctg_table.num_ctg_occs(contig_id)
    }

    fn get_base_index(&self) -> &BaseIndex {
        &self.base
    }

    delegate! {
        to self.ctg_table {
            fn get_ref_name(&self, id: usize) -> Option<&String> ;
            fn get_ref_names(&self) -> &[String];
            fn get_encoded_contig_occs(&self, contig_id: usize) -> &[u64];
        }
    }

    fn to_info(&self) -> Info {
        let base = &self.base;
        Info {
            sampling_type: PufferfishType::DenseRS,
            index_version: base.index_version,
            reference_gfa: base.reference_gfa.clone(),
            kmer_size: base.k,
            num_kmers: base.num_kmers(),
            num_contigs: base.num_contigs(),
            seq_len: base.seq_len(),
            have_ref_seq: base.ref_seq.is_some(),
            have_edge_vec: base.have_edge_vec,
            seq_hash: base.seq_hash.clone(),
            name_hash: base.name_hash.clone(),
            seq_hash_512: base.seq_hash_512.clone(),
            name_hash_512: base.name_hash_512.clone(),
            decoy_seq_hash: base.decoy_seq_hash.clone(),
            decoy_name_hash: base.decoy_name_hash.clone(),
            num_decoys: base.num_decoys,
            sample_size: None, // Putting into enum would break c++ compatiblity)
            extension_size: None,
            first_decoy_index: base.first_decoy_index,
            keep_duplicates: base.keep_duplicates,
        }
    }

    fn iter_refs(&self) -> RefIterator<Self> {
        // fn iter_refs<'a: 'b, 'b>(&'a self) -> RefIterator<'b, Self> {
        if !self.has_ref_seq() {
            // TODO: We may want to move this panic inside of the kmer iterators
            // if the BaseIndexReference supports more functionalities that don't
            // require the sequence
            panic!("Cannot iterate over references with no stored refseq")
        }
        RefIterator::new(self)
    }
}

impl<'a> PuffQuery<'a> for DenseIndex {
    type HitsT = ProjectedHits<'a>;

    fn get_ref_pos(&'a self, kmer: &CanonicalKmer) -> Option<Self::HitsT> {
        self.check_query(&kmer.get_fw_mer());
        // `?` operator is ok on Option that impls Try
        // returns none if none
        // position of match on seq
        let (pos, mapped_orientation) = self.get_kmer_pos(kmer)?;
        let pos = pos as usize;
        let contig_id = self.base.bv.rank(pos) as usize;

        let ctg_start_pos = self.contig_start_pos(contig_id);
        let contig_len = self.contig_len(contig_id);
        let offset = pos - ctg_start_pos;

        let hits = self.get_encoded_contig_occs(contig_id);
        Some(ProjectedHits {
            kmer: kmer.get_fw_mer(),
            mapped_orientation,
            offset,
            hits,
            contig_len,
            contig_id,
        })
    }
}

impl<'a> CachedQuery<'a, QueryCache> for DenseIndex {
    type HitsT = ProjectedHits<'a>;
    fn get_ref_pos_with_cache(
        &'a self,
        kmer: &CanonicalKmer,
        qc: &mut QueryCache,
    ) -> Option<Self::HitsT> {
        self.check_query(&kmer.get_fw_mer());
        // position of match on seq
        let (pos, mapped_orientation) = self.get_kmer_pos(kmer)?;

        let same_contig_as_cache =
            (qc.contig_start <= pos) && (pos < qc.contig_start + qc.contig_len);

        let contig_id;
        let offset;
        let contig_len;

        if same_contig_as_cache {
            contig_id = qc.contig_id;
            contig_len = qc.contig_len;
            offset = pos - qc.contig_start;
        } else {
            contig_id = self.base.bv.rank(pos) as usize;

            // start position of contig on seq
            let ctg_start_pos = self.contig_start_pos(contig_id);
            contig_len = self.contig_len(contig_id);
            offset = pos - ctg_start_pos;

            qc.contig_len = contig_len;
            qc.contig_id = contig_id;
            qc.contig_start = ctg_start_pos;
        }

        let hits = self.get_encoded_contig_occs(contig_id);
        Some(ProjectedHits {
            kmer: kmer.get_fw_mer(),
            mapped_orientation,
            offset,
            hits,
            contig_len,
            contig_id,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::PufferfishType;
    use crate::test_utils::*;

    #[test]
    fn pufferfish_type() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        assert_eq!(pi.to_info().sampling_type, PufferfishType::DenseRS);
    }
    #[test]
    fn tiny_contig_offsets() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        assert_eq!(pi.num_contigs(), pi.ctg_table.contig_offsets.len() - 1);
    }

    #[test]
    fn yeast_contig_offsets() {
        let p = to_abs_path(YEAST_CHR01_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        assert_eq!(pi.num_contigs(), pi.ctg_table.contig_offsets.len() - 1);
    }

    #[test]
    fn tiny_query_does_not_exist() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        let kmers = vec!["tat", "ata", "act", "ctg", "cct"];
        let kmers = kmers.iter().map(|s| CanonicalKmer::from(*s));

        for k in kmers {
            assert_eq!(pi.get_kmer_pos(&k), None);
        }
    }

    /******************************************************************************/
    // Tiny example: invalid queries panic
    /******************************************************************************/
    #[test]
    #[should_panic]
    fn tiny_kmer_too_big() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        pi.get_kmer_pos(&CanonicalKmer::from("aaaaaa"));
    }

    #[test]
    #[should_panic]
    fn tiny_kmer_too_small() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
        pi.get_kmer_pos(&CanonicalKmer::from("aa"));
    }

    /******************************************************************************/
    // Tiny example: check kmer positions on seq
    /******************************************************************************/
    #[test]
    fn test_tiny_kmer_positions() {
        let p = to_abs_path(TINY_INDEX);
        let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

        // Indexed string: // AAACCC
        let aaa_pos = 0;
        let aac_pos = 1;
        let acc_pos = 2;
        let ccc_pos = 3;

        // binary representations;
        let mut aaa = CanonicalKmer::from("aaa");
        let mut aac = CanonicalKmer::from("aac");
        let mut acc = CanonicalKmer::from("acc");
        let mut ccc = CanonicalKmer::from("ccc");

        assert_eq!(
            pi.get_kmer_pos(&aaa),
            Some((aaa_pos, MappedOrientation::Forward))
        );
        assert_eq!(
            pi.get_kmer_pos(&aac),
            Some((aac_pos, MappedOrientation::Forward))
        );
        assert_eq!(
            pi.get_kmer_pos(&acc),
            Some((acc_pos, MappedOrientation::Forward))
        );
        assert_eq!(
            pi.get_kmer_pos(&ccc),
            Some((ccc_pos, MappedOrientation::Forward))
        );

        aaa.swap();
        assert_eq!(
            pi.get_kmer_pos(&aaa),
            Some((aaa_pos, MappedOrientation::Backward))
        );
        aac.swap();
        assert_eq!(
            pi.get_kmer_pos(&aac),
            Some((aac_pos, MappedOrientation::Backward))
        );
        acc.swap();
        assert_eq!(
            pi.get_kmer_pos(&acc),
            Some((acc_pos, MappedOrientation::Backward))
        );
        ccc.swap();
        assert_eq!(
            pi.get_kmer_pos(&ccc),
            Some((ccc_pos, MappedOrientation::Backward))
        );
    }
}
