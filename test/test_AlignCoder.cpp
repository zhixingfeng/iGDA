//
//  test_AlignCoder.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/misc.h"
#include "../src/modules/modules.h"
TEST_CASE("test AlignCoderSNV::encode", "[hide]")
{
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->encode("../data/MSSA_61_forward.m5", "../results/MSSA_61_forward_encode_snv.txt");
}

TEST_CASE("test AlignCoderSNV::encode negative strand","[hide]")
{
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->encode("../data/B_10_cons.m5", "../results/B_10_cons.encode");
}


TEST_CASE("test AlignCoderSNV::decode", "[hide]")
{
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;

    REQUIRE(p_aligncoder->decode(2).first==1); REQUIRE(p_aligncoder->decode(2).second=='C');
    REQUIRE(p_aligncoder->decode(3).first==1); REQUIRE(p_aligncoder->decode(3).second=='G');
    REQUIRE(p_aligncoder->decode(4).first==1); REQUIRE(p_aligncoder->decode(4).second=='T');
    REQUIRE(p_aligncoder->decode(5).first==2); REQUIRE(p_aligncoder->decode(5).second=='A');
    REQUIRE(p_aligncoder->decode(6).first==2); REQUIRE(p_aligncoder->decode(6).second=='C');
    REQUIRE(p_aligncoder->decode(7).first==2); REQUIRE(p_aligncoder->decode(7).second=='G');
    REQUIRE(p_aligncoder->decode(8).first==2); REQUIRE(p_aligncoder->decode(8).second=='T');
    
    REQUIRE(p_aligncoder->decode(3894).first==974); REQUIRE(p_aligncoder->decode(3894).second=='C');
    REQUIRE(p_aligncoder->decode(372).first==93);REQUIRE(p_aligncoder->decode(372).second=='T');
    
    
    
}





TEST_CASE("test AlignCoderSNV::encode sam format", "[hide]")
{
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderSam alignreadersam;
    alignreadersam.getref("../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.ref.fa");
    p_aligncoder->setAlignReader(&alignreadersam);
    p_aligncoder->encode("../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.sam",
                         "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.encode.fromsam");
    
    string cmd = "sort ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.encode > ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.encode.sorted";
    system(cmd.c_str());
    
    cmd = "sort ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.encode.fromsam > ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.encode.fromsam.sorted";
    system(cmd.c_str());
    
}

TEST_CASE("test AlignCoderSNV::recode()","[hide]")
{
    //string m5_file = "../results/recode/m170701_065929_00127_c101204132550000001823285910241766_s1_p0.m5";
    //string var_file = "../results/recode/m170701_065929_00127_c101204132550000001823285910241766_s1_p0.var";
    //string recode_file = "../results/recode/m170701_065929_00127_c101204132550000001823285910241766_s1_p0.recode";
    
    string m5_file = "../results/pt_recode/align_to_consensus_trim.100.m5";
    string var_file = "../results/pt_recode/align_to_consensus_trim.var";
    string recode_file = "../results/pt_recode/align_to_consensus_trim.100.recode";

    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->recode(m5_file, var_file, recode_file, 10, 10, true);
    
    
}



TEST_CASE("test AlignCoderSNV::recode() debug", "[hide]")
{
    string m5_file = "../results/recode_debug/ERR2672423_tig00000206_quiver_line_3842.m5";
    string var_file = "../results/recode_debug/ERR2672423_tig00000206_quiver.var";
    string recode_file = "../results/recode_debug/ERR2672423_tig00000206_quiver.recode";
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->recode(m5_file, var_file, recode_file, 10, 10, true);
    
    
}

TEST_CASE("test AlignCoderSNV::recode() check putative errors", "[hide]")
{
    string m5_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.m5";
    string var_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.var";
    string recode_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.recode.debug";
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->recode(m5_file, var_file, recode_file, 10, 10, true);
    //p_aligncoder->recode_legacy(m5_file, var_file, recode_file, 10, 10, true);
}


TEST_CASE("test AlignCoderSNV::encode() encode reference", "[hide]")
{
    string m5_file = "../data/test_encode/realign.m5.1000";
    string encode_file = "../data/test_encode/realign.encode.1000";
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->encode(m5_file, encode_file);
}


TEST_CASE("test recode using masked reads", "[hide]")
{
    //string recode_file = "../data/nanopore_kp/igda_pipe_maskqv/merged_maskqv.debug.recode";
    string recode_file = "../data/nanopore_kp/igda_pipe_maskqv/merged.recode";
    string m5_file = "../data/nanopore_kp/igda_pipe_maskqv/merged_maskqv.m5";
    string var_file = "../data/nanopore_kp/igda_pipe_maskqv/realign_multi_nm.var";
    
    AlignCoderSNV aligncodersnv;
    AlignReaderM5 alignreaderm5;
    aligncodersnv.setAlignReader(&alignreaderm5);
    aligncodersnv.recode(m5_file, var_file, recode_file, 10, 10);

}

TEST_CASE("test multithrad recode()", "[hide]")
{
    int64_t nthread = 4;
    string align_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_recode_multithread/realign.m5";
    string var_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_recode_multithread/realign.var";
    string recode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_recode_multithread/realign.recode.multithread";
    AlignCoderSNV aligncodersnv;
    AlignReaderM5 alignreaderm5;
    AlignCoder *p_aligncoder = &aligncodersnv;
    
    //p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->recode_multithread(align_file, var_file, recode_file, 10, 10, 4);
    
    
    // split aligned reads
    /*vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    int64_t nlines = ceil ( (double)reads_range.size() / nthread);
    string cmd = "split -a 6 -l " + to_string(nlines) + " " + align_file + " " + align_file + ".tmp.split.part.";
    cout << cmd << endl; system(cmd.c_str());*/
    
    
    //split_file(align_file, align_file + ".c.tmp.split", nlines);
    
    /*// multithread recode
    for (auto i = 0; i < nthread; ++i){
        string cur_m5file = align_file + ".c.tmp.split.part." + to_string(i);
        string cur_outfile = align_file + ".c.tmp.split.part." + to_string(i) + ".recode";
        cmd = "igda recode " + cur_m5file + " " + var_file + " " + cur_outfile + " &";
        cout << cmd << endl; system(cmd.c_str());
        //system("igda recode " + )
    }
    system("wait");*/
    
}
TEST_CASE("test split_file", "[hide]")
{
    int64_t nthread = 4;
    string align_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_recode_multithread/realign.m5";
    cout << get_file_nlines(align_file) << endl;
    
}



