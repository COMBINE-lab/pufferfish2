use serde_json;
use simple_sds::ops::BitVec;
use std::fs::File;
use std::io::BufReader;

use super::contig::*;
use super::info::*;
use super::*;
use crate::boophf::MPHF;
use crate::cpp::DeserializeFromCpp;
use crate::test_utils::*;

#[test]
fn consts() {
    use consts::fp::*;
    let pf = PufferfishFilePaths::new("");
    assert_path_eq(pf.prefix, "");
    assert_path_eq(pf.complete_ref_lens, COMPLETE_REF_LENS);
    assert_path_eq(pf.ctable, CTABLE);
    assert_path_eq(pf.ctg_offsets, CTG_OFFSETS);
    assert_path_eq(pf.mphf, MPHF);
    assert_path_eq(pf.pos, POS);
    assert_path_eq(pf.rank, RANK);
    assert_path_eq(pf.ref_accum_lens, REF_ACCUM_LENS);
    assert_path_eq(pf.ref_lens, REF_LENS);
    assert_path_eq(pf.ref_seq, REF_SEQ);
    assert_path_eq(pf.seq, SEQ);
    assert_path_eq(pf.duplicate_clusters_tsv, DUPLICATE_CLUSTERS_TSV);
    assert_path_eq(pf.info_json, INFO_JSON);
    assert_path_eq(pf.ref_indexing_log, REF_INDEXING_LOG);
}

/******************************************************************************/
// Info
/******************************************************************************/
#[test]
fn info() {
    let p = to_abs_path(YEAST_CHR01_INDEX);
    let fps = PufferfishFilePaths::new(&p);
    let info_file = File::open(fps.info_json).unwrap();
    let rdr = BufReader::new(info_file);
    let info: Info = serde_json::from_reader(rdr).unwrap();

    let info_val = Info {
        index_version: 4,
        reference_gfa: vec![String::from("yeastchr01/pufferfish_index")],
        sampling_type: PufferfishType::Dense,
        kmer_size: 31,
        num_kmers: 221918,
        num_contigs: 577,
        seq_len: 239228,
        have_ref_seq: true,
        have_edge_vec: false,
        seq_hash: String::from("3c5c06b2ccb802798265a543cc6511d954a0a64a522c3f6af05be0553d6f0a62"),
        name_hash: String::from("7fcb45af27576afe48ccd78ec97f1f259a20ba5af218326e5ca0435002845803"),
        seq_hash_512: String::from("959cb1883fc1ca9ae1394ceb475a356ead1ecceff5824ae730c25c3c8cb58c90d823ca31107c6163cd4b32b033a1c1af84d11adc9271dd181d331404286184cb"),
        name_hash_512: String::from("78631229740885755bb5f7074b822b25fd339f473a34bddd7aab1ca67f2f2fbbd3501c9d2993a8c3bb0e83dc358943b359043e3a316b83ae5afe821f59c7ef66"),
        decoy_seq_hash: String::from("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
        decoy_name_hash: String::from("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
        extension_size: None,
        sample_size: None,
        num_decoys: 0,
        first_decoy_index: 18446744073709551615,
        keep_duplicates: false
    };
    assert_eq!(info, info_val);
}

#[test]
fn test_loaded_info_fields() {
    // TODO: maybe we can remove this test
    let p = to_abs_path(YEAST_CHR01_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    assert_eq!(pi.k(), 31);
    assert_eq!(pi.num_kmers(), 221918);
    assert_eq!(pi.num_contigs(), 577);
    assert_eq!(pi.seq_len(), 239228);

    assert_eq!(pi.base.num_decoys, 0);
    assert_eq!(pi.base.first_decoy_index, 18446744073709551615);
}

#[test]
fn tiny_to_info() {
    let p = to_abs_path(TINY_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);

    let mut info = Info::load(files.info_json);
    // change the loaded info to rust specific index type
    let loaded_info = pi.to_info();
    info.sampling_type = PufferfishType::DenseRS;

    assert_eq!(loaded_info, info);
}

#[test]
fn yeast_to_info() {
    let p = to_abs_path(YEAST_CHR01_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);

    let mut info = Info::load(files.info_json);
    // change the loaded info to rust specific index type
    // change the loaded info to rust specific index type
    let loaded_info = pi.to_info();
    info.sampling_type = PufferfishType::DenseRS;

    assert_eq!(loaded_info, info);
}

#[test]
fn test_tiny_getters() {
    let p = to_abs_path(TINY_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);
    let info = Info::load(files.info_json);
    assert_eq!(pi.num_kmers(), info.num_kmers);
    assert_eq!(pi.seq_len(), info.seq_len);
}

#[test]
fn test_yeast_getters() {
    let p = to_abs_path(TINY_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);
    let info = Info::load(files.info_json);
    assert_eq!(pi.num_kmers(), info.num_kmers);
    assert_eq!(pi.seq_len(), info.seq_len);
}

/******************************************************************************/
// Boundary bitvec
/******************************************************************************/
#[test]
fn tiny_bv_ones_equals_num_contigs() {
    let p = to_abs_path(TINY_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);
    let info = Info::load(files.info_json);

    assert_eq!(pi.base.bv.count_ones(), info.num_contigs);
    assert_eq!(pi.num_contigs(), info.num_contigs);
}

#[test]
fn yeast_bv_ones_equals_num_contigs() {
    let p = to_abs_path(YEAST_CHR01_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);
    let info = Info::load(files.info_json);

    assert_eq!(pi.base.bv.count_ones(), info.num_contigs);
    assert_eq!(pi.num_contigs(), info.num_contigs);
}

/******************************************************************************/
// Tiny example: check kmer positions on seq, querying the bits
//   - May be redundant
/******************************************************************************/
#[test]
fn test_tiny_u64_positions() {
    // Check pos.bin of "tiny" example that indexes the sequence:
    // AAACCC

    let p = to_abs_path(TINY_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();

    let aaa = 0;
    let pos_i = pi.base.mphf.lookup(&aaa).unwrap() as usize;
    let pos = pi.pos.get(pos_i);
    assert_eq!(pos, 0);

    // ordered from low-to-high order bits <=> left-to-right in cannonical kmer
    let aaa = 0b000000;
    let aac = 0b010000;
    let acc = 0b010100;
    let ccc = 0b010101;

    let aaa_i = pi.base.mphf.lookup(&aaa).unwrap() as usize;
    let aaa_pos = pi.pos.get(aaa_i);

    let aac_i = pi.base.mphf.lookup(&aac).unwrap() as usize;
    let aac_pos = pi.pos.get(aac_i);

    let acc_i = pi.base.mphf.lookup(&acc).unwrap() as usize;
    let acc_pos = pi.pos.get(acc_i);

    let ccc_i = pi.base.mphf.lookup(&ccc).unwrap() as usize;
    let ccc_pos = pi.pos.get(ccc_i);

    let aaa_pos_checked = 0;
    let aac_pos_checked = 1;
    let acc_pos_checked = 2;
    let ccc_pos_checked = 3;

    assert_eq!(aaa_pos, aaa_pos_checked);
    assert_eq!(aac_pos, aac_pos_checked);
    assert_eq!(acc_pos, acc_pos_checked);
    assert_eq!(ccc_pos, ccc_pos_checked);
}

/******************************************************************************/
// Tiny example: get_ref_pos -- full pufferfish query
/******************************************************************************/

#[test]
fn tiny_get_ref_pos() {
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

    let contig_info = ContigOcc {
        ref_id: 0,
        pos: 0,
        orientation: ContigOrientation::Forward,
    };
    let aaa_result = pi.get_ref_pos(&aaa).unwrap();
    assert_eq!(aaa_result.contig_id, 0);
    assert_eq!(aaa_result.contig_len, 6);
    assert_eq!(aaa_result.kmer, aaa.get_fw_mer());
    assert_eq!(aaa_result.len(), 1);
    let aaa_contig_info = ContigOcc::from(aaa_result.hits[0]);
    assert_eq!(aaa_contig_info, contig_info);
    assert_eq!(aaa_result.offset, aaa_pos);

    aaa.swap();
    let aaa_rc_result = pi.get_ref_pos(&aaa).unwrap();
    assert_eq!(aaa_rc_result.offset, aaa_result.offset);
    assert_eq!(aaa_rc_result.kmer, aaa_result.kmer.to_reverse_complement());
    assert_eq!(aaa_rc_result.hits, aaa_result.hits);

    let aac_result = pi.get_ref_pos(&aac).unwrap();
    assert_eq!(aac_result.contig_id, 0);
    assert_eq!(aac_result.contig_len, 6);
    assert_eq!(aac_result.kmer, aac.get_fw_mer());
    assert_eq!(aac_result.len(), 1);
    let aac_contig_info = ContigOcc::from(aac_result.hits[0]);
    assert_eq!(aac_contig_info, contig_info);
    assert_eq!(aac_result.offset, aac_pos);

    aac.swap();
    let aac_rc_result = pi.get_ref_pos(&aac).unwrap();
    assert_eq!(aac_rc_result.offset, aac_result.offset);
    assert_eq!(aac_rc_result.kmer, aac_result.kmer.to_reverse_complement());
    assert_eq!(aac_rc_result.hits, aac_result.hits);

    let acc_result = pi.get_ref_pos(&acc).unwrap();
    assert_eq!(acc_result.contig_id, 0);
    assert_eq!(acc_result.contig_len, 6);
    assert_eq!(acc_result.kmer, acc.get_fw_mer());
    assert_eq!(acc_result.len(), 1);
    let acc_contig_info = ContigOcc::from(acc_result.hits[0]);
    assert_eq!(acc_contig_info, contig_info);
    assert_eq!(acc_result.offset, acc_pos);

    acc.swap();
    let acc_rc_result = pi.get_ref_pos(&acc).unwrap();
    assert_eq!(acc_rc_result.offset, acc_result.offset);
    assert_eq!(acc_rc_result.kmer, acc_result.kmer.to_reverse_complement());
    assert_eq!(acc_rc_result.hits, acc_result.hits);

    let ccc_result = pi.get_ref_pos(&ccc).unwrap();
    assert_eq!(ccc_result.contig_id, 0);
    assert_eq!(ccc_result.contig_len, 6);
    assert_eq!(ccc_result.kmer, ccc.get_fw_mer());
    assert_eq!(ccc_result.len(), 1);
    let ccc_contig_info = ContigOcc::from(ccc_result.hits[0]);
    assert_eq!(ccc_contig_info, contig_info);
    assert_eq!(ccc_result.offset, ccc_pos);

    ccc.swap();
    let ccc_rc_result = pi.get_ref_pos(&ccc).unwrap();
    assert_eq!(ccc_rc_result.offset, ccc_result.offset);
    assert_eq!(ccc_rc_result.kmer, ccc_result.kmer.to_reverse_complement());
    assert_eq!(ccc_rc_result.hits, ccc_result.hits);
}

#[test]
fn tiny_get_ref_pos_with_cache() {
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

    let contig_info = ContigOcc {
        ref_id: 0,
        pos: 0,
        orientation: ContigOrientation::Forward,
    };

    let mut qc = QueryCache {
        contig_start: usize::MAX,
        contig_len: usize::MAX,
        contig_id: usize::MAX,
    };

    let aaa_result = pi.get_ref_pos_with_cache(&aaa, &mut qc).unwrap();
    assert_eq!(aaa_result.contig_id, 0);
    assert_eq!(aaa_result.contig_len, 6);
    assert_eq!(aaa_result.kmer, aaa.get_fw_mer());
    assert_eq!(aaa_result.len(), 1);
    let aaa_contig_info = ContigOcc::from(aaa_result.hits[0]);
    assert_eq!(aaa_contig_info, contig_info);
    assert_eq!(aaa_result.offset, aaa_pos);

    aaa.swap();
    let aaa_rc_result = pi.get_ref_pos_with_cache(&aaa, &mut qc).unwrap();
    assert_eq!(aaa_rc_result.offset, aaa_result.offset);
    assert_eq!(aaa_rc_result.kmer, aaa_result.kmer.to_reverse_complement());
    assert_eq!(aaa_rc_result.hits, aaa_result.hits);

    let aac_result = pi.get_ref_pos_with_cache(&aac, &mut qc).unwrap();
    assert_eq!(aac_result.contig_id, 0);
    assert_eq!(aac_result.contig_len, 6);
    assert_eq!(aac_result.kmer, aac.get_fw_mer());
    assert_eq!(aac_result.len(), 1);

    aac.swap();
    let aac_contig_info = ContigOcc::from(aac_result.hits[0]);
    assert_eq!(aac_contig_info, contig_info);
    assert_eq!(aac_result.offset, aac_pos);
    let aac_rc_result = pi.get_ref_pos_with_cache(&aac, &mut qc).unwrap();
    assert_eq!(aac_rc_result.offset, aac_result.offset);
    assert_eq!(aac_rc_result.kmer, aac_result.kmer.to_reverse_complement());
    assert_eq!(aac_rc_result.hits, aac_result.hits);

    let acc_result = pi.get_ref_pos_with_cache(&acc, &mut qc).unwrap();
    assert_eq!(acc_result.contig_id, 0);
    assert_eq!(acc_result.contig_len, 6);
    assert_eq!(acc_result.kmer, acc.get_fw_mer());
    assert_eq!(acc_result.len(), 1);

    acc.swap();
    let acc_contig_info = ContigOcc::from(acc_result.hits[0]);
    assert_eq!(acc_contig_info, contig_info);
    assert_eq!(acc_result.offset, acc_pos);
    let acc_rc_result = pi.get_ref_pos_with_cache(&acc, &mut qc).unwrap();
    assert_eq!(acc_rc_result.offset, acc_result.offset);
    assert_eq!(acc_rc_result.kmer, acc_result.kmer.to_reverse_complement());
    assert_eq!(acc_rc_result.hits, acc_result.hits);

    let ccc_result = pi.get_ref_pos_with_cache(&ccc, &mut qc).unwrap();
    assert_eq!(ccc_result.contig_id, 0);
    assert_eq!(ccc_result.contig_len, 6);
    assert_eq!(ccc_result.kmer, ccc.get_fw_mer());
    assert_eq!(ccc_result.len(), 1);

    ccc.swap();
    let ccc_contig_info = ContigOcc::from(ccc_result.hits[0]);
    assert_eq!(ccc_contig_info, contig_info);
    assert_eq!(ccc_result.offset, ccc_pos);
    let ccc_rc_result = pi.get_ref_pos_with_cache(&ccc, &mut qc).unwrap();
    assert_eq!(ccc_rc_result.offset, ccc_result.offset);
    assert_eq!(ccc_rc_result.kmer, ccc_result.kmer.to_reverse_complement());
    assert_eq!(ccc_rc_result.hits, ccc_result.hits);
}
/******************************************************************************/
// Contig table loading
/******************************************************************************/
#[test]
fn tiny_contig_table() {
    let p = to_abs_path(TINY_INDEX);
    let pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);

    let (names, exts, ctable) =
        DenseContigTable::deserialize_from_cpp_helper(&files.ctable).unwrap();
    assert_eq!(names, vec!["tiny"]);
    assert_eq!(exts.len(), 1);

    assert_eq!(ctable.len(), pi.num_contigs());
    assert_eq!(ctable, pi.ctg_table.ctable);

    assert_eq!(ctable.len(), 1);
    let info = ContigOcc::from(ctable[0]);
    let info_correct = ContigOcc {
        pos: 0,
        ref_id: 0,
        orientation: ContigOrientation::Forward,
    };
    assert_eq!(info, info_correct);
}

#[test]
fn yeast_contig_table() {
    // TODO check more than just the length of names and extensions
    // The only "good" way to do this is to run the internal validate function
    let p = to_abs_path(YEAST_CHR01_INDEX);
    let _pi = DenseIndex::deserialize_from_cpp(&p).unwrap();
    let files = PufferfishFilePaths::new(&p);
    let _info = Info::load(files.info_json);

    let (names, exts, _ctable) =
        DenseContigTable::deserialize_from_cpp_helper(&files.ctable).unwrap();
    assert_eq!(names.len(), 1);
    assert_eq!(exts.len(), 1);
}
