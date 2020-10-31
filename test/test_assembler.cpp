//
//  test_assemble.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/io.h"
#include "../src/misc/basic.h"
#include "../src/modules/assemble/assembler.h"
#include "../src/modules/dforest/dforestsnvstl.h"


TEST_CASE("test assembler::get_variants()", "[hide]")
{
    string dforest_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.dforest.uniq.max";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    Assembler assembler;
    assembler.get_variants(dforest_file, out_file, 0.8);
    
}

TEST_CASE("test assembler::reduce_dim()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode";
    string var_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim";
    Assembler assembler;
    assembler.reduce_dim(encode_file, var_file, out_file);
}

TEST_CASE("test assembler::dist()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.dist";
    Assembler assembler;
    assembler.dist(encode_file, m5_file, out_file);
}

TEST_CASE("test assembler::dist_rdim() legacy", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string var_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.dist.withvar";
    Assembler assembler;
    assembler.dist_rdim(encode_file, m5_file, var_file, out_file);
}

TEST_CASE("test assembler::jaccard_index()", "[hide]")
{
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.jaccard.5000";
    Assembler assembler;
    assembler.jaccard_index(encode_file, m5_file, out_file, 0.5);
}

TEST_CASE("test assembler::jaccard_index_min()", "[hide]")
{
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.jaccard_min.5000";
    Assembler assembler;
    assembler.jaccard_index_min(encode_file, m5_file, out_file, 0.1);
}


TEST_CASE("test assembler::haplo_seq_construct()", "[hide]")
{
    string encode_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var.encode";
    string ref_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.ref";
    string out_dir = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642_haplotype";
    
    // read encodefile
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // read reference file
    unordered_map<string, string> ref_seq_all;
    read_fasta(ref_file, ref_seq_all);
    if (ref_seq_all.size() != 1)
        throw runtime_error("number of chrosome is not 1.");
    string ref_seq = ref_seq_all.begin()->second;
    

    // construct haplotype
    for (int i = 0; i < (int) encode_data.size(); ++i){
        Assembler assembler;
        string haplo_seq;
        assembler.haplo_seq_construct(encode_data[i], ref_seq, haplo_seq);
        
        string outfile = out_dir + "/haplotype_" + to_string(i) + ".fa";
        ofstream fs_outfile;
        open_outfile(fs_outfile, outfile);
        fs_outfile << ">haplotype_" + to_string(i) << endl;
        fs_outfile << haplo_seq << endl;
        fs_outfile.close();
    }
    
}

TEST_CASE("test assembler::dist_rdim()", "[hide]")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.5000";
    string m5_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    string var_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    string out_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.dist.withvar";
    Assembler assembler;
    assembler.dist_rdim(encode_file, m5_file, var_file, out_file);
}


TEST_CASE("test assembler::ann_clust", "[hide]")
{
    //string recode_file = "../results/pt_recode/align_to_consensus_trim.5000.recode";
    //string recode_ref_file = "../results/pt_recode/align_to_consensus_trim.5000.recode.ref";
    //string m5_file = "../results/pt_recode/align_to_consensus_trim.5000.m5";
    //string var_file = "../results/pt_recode/align_to_consensus_trim.var";
    //string out_file = "../results/pt_recode/align_to_consensus_trim.ann";
    
    string encode_file = "../results/pt_rpoBC/igda/realign.encode.rdim";
    string recode_file = "../results/pt_rpoBC/igda/realign.recode";
    string recode_ref_file = "../results/pt_rpoBC/igda/realign.recode.ref";
    string m5_file = "../results/pt_rpoBC/igda/realign.m5";
    string var_file = "../results/pt_rpoBC/igda/realign.var";
    string out_file = "../results/pt_rpoBC/igda/realign.ann.debug";
    
    Assembler assembler;
    
    assembler.ann_clust_recode(recode_file, recode_ref_file, encode_file, m5_file, var_file, 10, 0.2, 0.8, 23, 46, 0.5);
    
    //assembler.ann_clust_recode_legacy(recode_file, recode_ref_file, m5_file, var_file, 20, 0.2, 0.7, 20, 50, 0.02);
    //assembler.ann_clust(encode_file, m5_file, var_file, 20, 0.2, 0.7, 30, 200, 0.02);
    
    
    vector<int64_t> idx;
    assembler.find_nccontigs(idx);
    assembler.print_rl_ann_clust(out_file + ".igda_tmp", true, idx);
    string cmd = "sort -u -s -k2,2n -k3,3n -k1,1 " + out_file + ".igda_tmp" + " > " + out_file;
    cout << cmd << endl; system(cmd.c_str());
    //cout << idx << endl;
    
}

TEST_CASE("test assembler::ann_clust, debug pt_recode", "[hide]")
{
    string recode_file = "../results/pt_ann_recode_debug/align_to_consensus_trim.recode";
    string recode_ref_file = "../results/pt_ann_recode_debug/align_to_consensus_trim.recode.ref";
    string m5_file = "../results/pt_ann_recode_debug/align_to_consensus_trim.m5";
    string var_file = "../results/pt_ann_recode_debug/align_to_consensus_trim.var";
    string out_file = "../results/pt_ann_recode_debug/align_to_consensus_trim.ann";
    
    Assembler assembler;
    //assembler.ann_clust_recode_legacy(recode_file, recode_ref_file, m5_file, var_file, 12, 0.2, 0.8, 20, 50, 0.5);
    //assembler.ann_clust(encode_file, m5_file, var_file, 20, 0.2, 0.7, 30, 200, 0.02);
    
    
    vector<int64_t> idx;
    assembler.find_nccontigs(idx);
    assembler.print_rl_ann_clust(out_file, false, idx);
    cout << idx << endl;
    
    //assembler.print_rl_ann_clust(out_file+".seq", true);
}


TEST_CASE("test assembler::ann_clust compare new and legacy recoding algorithm", "[hide]")
{
    string recode_file = "../results/pt_ann_recode_cmp_new_and_legacy/align_to_consensus_trim.recode";
    string recode_ref_file = "../results/pt_ann_recode_cmp_new_and_legacy/align_to_consensus_trim.recode.ref";
    string m5_file = "../results/pt_ann_recode_cmp_new_and_legacy/align_to_consensus_trim.m5";
    string var_file = "../results/pt_ann_recode_cmp_new_and_legacy/align_to_consensus_trim.var";
    string out_file = "../results/pt_ann_recode_cmp_new_and_legacy/align_to_consensus_trim.ann";
    
    Assembler assembler;
    //assembler.ann_clust_recode_legacy(recode_file, recode_ref_file, m5_file, var_file, 12, 0.2, 0.8, 20, 50, 0.5);

}

TEST_CASE("test assembler::assign_reads_to_contigs()", "[hide]")
{
    string ann_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.ann.tested.ft";
    string recode_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.recode";
    string m5_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.m5";
    //string ann_file = "../results/pt_rpoBC/igda/realign.ann.tested.ft";
    //string recode_file = "../results/pt_rpoBC/igda/realign.recode";
    //string m5_file = "../results/pt_rpoBC/igda/realign.m5";
    
    Assembler assembler;
    cout << "load ann" << endl;
    assembler.read_ann_results(ann_file);
    
    cout << "load recode_data" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    cout << "load reads_range" << endl;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, m5_file);
    
    assembler.assign_reads_to_contigs(recode_data, reads_range);
    
    assembler.print_rl_ann_clust(ann_file + ".count", true);
    
    Assembler assembler_2;
    assembler_2.read_ann_results(ann_file + ".count");
    assembler_2.print_rl_ann_clust(ann_file + ".count.copy", true);
    int tmp = -1;
    
}

TEST_CASE("test assembler::test_contigs()", "[hide]")
{
    string ref_file = "../results/pt_ann_assign_reads/consensus.fasta";
    string ann_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann";
    string recode_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode";
    string recode_ref_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ref";
    string m5_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.m5";
    //string out_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.test";
    
    Assembler assembler;
    cout << "load ref_file" << endl;
    assembler.load_homo_blocks(ref_file);
    
    cout << "load ann" << endl;
    assembler.read_ann_results(ann_file);
    
    cout << "load recode_data" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    cout << "load recode_ref_data" << endl;
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_ref_file);
    
    cout << "load reads_range" << endl;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, m5_file);
    
    
    assembler.test_contigs(recode_data, recode_ref_data, reads_range);
    assembler.print_rl_ann_clust("../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested", true);
    
    int x = 0;
}

TEST_CASE("test assembler::test_contigs_pairwise()", "[hide]")
{
    string ann_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.ann.tested.ft.count";
    string recode_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.recode";
    string ref_file = "../results/LoFreq_HBV/cat_wild_large.fasta";
    //string ann_file = "../results/pt_ecoli/igda/realign.ann.tested.ft.count";
    //string recode_file = "../results/pt_ecoli/igda/realign.recode";
    //string ref_file = "../results/pt_ecoli/ecoli_K12_MG1655.fasta";
    Assembler assembler;
    assembler.load_homo_blocks(ref_file);
    assembler.test_contigs_pairwise(ann_file, recode_file, ann_file + ".ft");
    
}

TEST_CASE("test assembler::read_ann_results()", "[hide]")
{
    string ann_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann";
    Assembler assembler;
    assembler.read_ann_results(ann_file);
    assembler.print_rl_ann_clust(ann_file + ".copy", true);
}


TEST_CASE("test assembler::find_nccontigs()", "[hide]")
{
    string ann_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann";
    Assembler assembler;
    assembler.read_ann_results(ann_file);
    vector<int64_t> idx;
    assembler.find_nccontigs(idx);
    cout << idx << endl;
}

TEST_CASE("test assembler::filter_ann()", "[hide]")
{
    Assembler assembler;
    assembler.filter_ann("../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested");
}


TEST_CASE("read_dot_file", "[hide]")
{
    cout << "original graph" << endl;
    Graph gp;
    read_dot_file(gp, "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft.dot");
    boost::print_graph(gp);
    
    cout << "transitive reducted graph" << endl;
    Graph gp_tred;
    read_dot_file(gp_tred, "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft.transitive_reduction.dot");
    boost::print_graph(gp_tred);
    
}

TEST_CASE("test Assembler::assemble", "[hide]")
{
    Assembler assembler;
    Graph gp;
    assembler.read_ann_results("../results/sa/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.ann.tested.ft");
    read_dot_file(gp, "../results/sa/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.ann.tested.ft.tred.dot");
    assembler.assemble(gp, "../results/sa/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.ann.tested.ft.assembled");
    
}


TEST_CASE("test assembler::test_contigs() (hbv)","[hide]")
{
    string ref_file = "../results/LoFreq_HBV/cat_wild_large.fasta";
    string ann_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.ann";
    string recode_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.recode";
    string recode_ref_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.recode.ref";
    string m5_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.m5";
    string out_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.ann.tested.debug";
    
    Assembler assembler;
    cout << "load ref_file" << endl;
    assembler.load_homo_blocks(ref_file);
    
    cout << "load ann" << endl;
    assembler.read_ann_results(ann_file);
    
    cout << "load recode_data" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    cout << "load recode_ref_data" << endl;
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_ref_file);
    
    cout << "load reads_range" << endl;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, m5_file);
    
    assembler.test_contigs(recode_data, recode_ref_data, reads_range);
    assembler.print_rl_ann_clust(out_file, true);
    
    int x = 0;
}

TEST_CASE("test assembler::ann_clust(), HBV", "[hide]")
{
    string recode_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.recode";
    string recode_ref_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.recode.ref";
    string encode_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.encode.rdim";
    string m5_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.m5";
    string var_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.var";
    string out_file = "../results/LoFreq_HBV/igda_result_large/m90min_1389bp_11l_aligned_reads/realign.ann.debug";
    
    Assembler assembler;
    
    assembler.ann_clust_recode(recode_file, recode_ref_file, encode_file, m5_file, var_file, 10, 0.2, 0.8, 15, 30, 0.5);
    
    vector<int64_t> idx;
    assembler.find_nccontigs(idx);
    assembler.print_rl_ann_clust(out_file + ".igda_tmp", true, idx);
    
    string cmd = "sort -u -s -k2,2n -k3,3n -k1,1 " + out_file + ".igda_tmp" + " > " + out_file;
    cout << cmd << endl; system(cmd.c_str());
    //cout << idx << endl;
    
    //assembler.print_rl_ann_clust(out_file+".seq", true);
}


TEST_CASE("test Assembler::assemble (HBV)", "[hide]")
{
    Assembler assembler;
    Graph gp;
    assembler.read_ann_results("../results/LoFreq_HBV/igda_result_large/3l/realign.ann.tested.ft");
    read_dot_file(gp, "../results/LoFreq_HBV/igda_result_large/3l/realign.ann.tested.ft.tred.dot");
    assembler.assemble(gp, "../results/LoFreq_HBV/igda_result_large/3l/realign.ann.tested.ft.assembled.debug");
    
}

TEST_CASE("test Assembler::ann_clust_recode (which reads correction)", "[hide]")
{
    string recode_file = "../data/test_correct_reads/SRR8054532.recode";
    string recode_ref_file = "../data/test_correct_reads/SRR8054532.recode.ref";
    string encode_file = "../data/test_correct_reads/SRR8054532.encode";
    string m5_file = "../data/test_correct_reads/SRR8054532.m5";
    string var_file = "../data/test_correct_reads/SRR8054532.var";
    string ann_file = "../data/test_correct_reads/SRR8054532.ann";
    
    Assembler assembler;
    assembler.ann_clust_recode(recode_file, recode_ref_file, encode_file, m5_file, var_file, 10, 0.3, 0.7, 25, 50, 0.5, true);
    vector<int64_t> idx;
    assembler.find_nccontigs(idx);
    assembler.print_rl_ann_clust(ann_file + ".igda_tmp", true, idx);
    string cmd = "sort -u -s -k2,2n -k3,3n -k1,1 " + ann_file + ".igda_tmp" + " > " + ann_file;
    cout << cmd << endl; system(cmd.c_str());
}


TEST_CASE("test assembler::test_contigs() (relative risk)", "[hide]")
{
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/k_pneumoniae_hs11286_chr.fasta";
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/debug_realign_multi_nm.ann";
    string recode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/realign_multi_nm.recode";
    string recode_ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/realign_multi_nm.recode.ref";
    string m5_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/realign.m5";
    string out_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/debug_realign_multi_nm.ann.tested";
    
    Assembler assembler;
    cout << "load ref_file" << endl;
    assembler.load_homo_blocks(ref_file);
    
    cout << "load ann" << endl;
    assembler.read_ann_results(ann_file);
    
    cout << "load recode_data" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    cout << "load recode_ref_data" << endl;
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_ref_file);
    
    cout << "load reads_range" << endl;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, m5_file);
    
    assembler.set_null_betadist(0.6466251, 21.90139);
    assembler.test_contigs(recode_data, recode_ref_data, reads_range);
    assembler.print_rl_ann_clust(out_file, true);
    
    assembler.filter_ann(out_file, 5, 10, 5);
    
    int x = 0;
}


TEST_CASE("test assembler::test_contigs_pairwise() (relative risk)", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/debug_realign_multi_nm.ann.tested.ft.count";
    string recode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/realign_multi_nm.recode";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/igda_pipe/k_pneumoniae_hs11286_chr.fasta";
    Assembler assembler;
    
    assembler.set_null_betadist(0.6466251, 21.90139);
    assembler.load_homo_blocks(ref_file);
    assembler.test_contigs_pairwise(ann_file, recode_file, ann_file + ".ft.debug", 5, 3, 10, 5);
    
}

TEST_CASE("test assembler::ann_clust_recode() compare", "[hide]")
{
    string recode_file = "../data/nanopore_kp/phase/c_1_t_15_m_30/realign.recode";
    string recode_ref_file = "../data/nanopore_kp/phase/c_1_t_15_m_30/realign.recode.ref";
    string encode_file = "../data/nanopore_kp/phase/c_1_t_15_m_30/realign.encode.rdim";
    string align_file = "../data/nanopore_kp/optimum_detect/realign.m5";
    string var_file = "../data/nanopore_kp/optimum_detect/realign.var";

    int mincvg = 1;
    double minprop = 0.2;
    double maxprop = 0.8;
    int topn = 15;
    int maxnn = 30;
    double minjaccard = 0.5;
    bool iscorrect = false;
    bool ishang = false;
    
    Assembler assembler;
    assembler.ann_clust_recode(recode_file, recode_ref_file, encode_file, align_file, var_file, mincvg, minprop, maxprop, topn, maxnn, minjaccard, iscorrect, ishang);
    assembler.print_rl_ann_clust("../data/nanopore_kp/phase/c_1_t_15_m_30/realign.ann.raw", true);
    
    mincvg = 10;
    assembler.ann_clust_recode(recode_file, recode_ref_file, encode_file, align_file, var_file, mincvg, minprop, maxprop, topn, maxnn, minjaccard, iscorrect, ishang);
    assembler.print_rl_ann_clust("../data/nanopore_kp/phase/c_10_t_15_m_30/realign.ann.raw", true);
    
    int x = 0;
    //assembler.ann_clust_recode(recodefileArg.getValue(), recodefileArg.getValue() + ".ref", encodefileArg.getValue(), alignfileArg.getValue(), varfileArg.getValue(), mincvgArg.getValue(),
    //                           minpropArg.getValue(), maxpropArg.getValue(), topnArg.getValue(), maxnnArg.getValue(), minjaccardArg.getValue(), iscorrectArg.getValue(), ishangArg.getValue());
}

TEST_CASE("itest assembler::correct_contigs()", "[hide]")
{
    string ann_file = "../data/nanopore_kp/phase/c_10_t_15_m_30/realign.ann.tested.ft.count.ft";
    string out_file = "../data/nanopore_kp/phase/c_10_t_15_m_30/realign.ann.tested.ft.count.ft.ct";
    Assembler assembler;
    assembler.correct_contigs(ann_file, out_file);
    
}

TEST_CASE("debug ann_to_graph", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_assemble/realign.ann.tested.ft.count.ft";
    Assembler assembler;
    Graph gp;
    assembler.ann_to_graph(gp, ann_file);
    ofstream fs_graph(ann_file + ".dot");
    boost::write_graphviz(fs_graph, gp);
    fs_graph.close();
}

TEST_CASE("test polish()", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/clpX_1/realign.ann.tested.ft.count.ft";
    string encode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/clpX_1/realign.encode.rdim";
    string m5_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/detect/clpX_1/realign.m5";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/Borrelia_burgdorferi_N40_chr.fna";
    string tmp_dir = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/detect/clpX_1/tmp_polish";
    string out_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/clpX_1/realign.ann.tested.ft.count.ft.polished";
    
    Assembler assembler;
    assembler.polish(ann_file, encode_file, m5_file, ref_file, out_file, tmp_dir, 0.3, 25, 15);

}
TEST_CASE("test dforest in polish()", "[hide]")
{
    //string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/uvrA_1/realign.ann.tested.ft.count.ft";
    string encode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/uvrA_1/realign.encode.rdim";
    string m5_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/detect/uvrA_1/realign.m5";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/Borrelia_burgdorferi_N40_chr.fna";
    string cmpreads_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/uvrA_1/realign.ann.tested.ft.count.ft.unpolished.cmpreads.test";
    string tmp_dir = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/detect/uvrA_1/tmp_polish";
    string out_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/uvrA_1/realign.ann.tested.ft.count.ft.dforest.test";
    
    AlignReaderM5 alignreaderm5;
    AlignCoderSNV aligncoder;
        
    DForestSNVSTL forestsnvstxxl(&alignreaderm5, &aligncoder);
    
    forestsnvstxxl.load_homo_blocks(ref_file);
    forestsnvstxxl.run(encode_file, m5_file, cmpreads_file, out_file, tmp_dir, 10, 1000, 1, 0.3, 1, 1, true);
    
}

TEST_CASE("debug Assemble::ann_to_graph()", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_polish/phase/clpX_1/realign.ann.tested.ft.count.ft";
    Assembler assembler;
    Graph gp;
    assembler.ann_to_graph(gp, ann_file);
    
    ofstream fs_graph(ann_file + ".dot");
    boost::write_graphviz(fs_graph, gp);
}

TEST_CASE("test ram issue of Assemble::assemble()", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_assemble_ram_issue/NC_001318.1_gene_661/realign.ann.tested.ft.count.ft";
    string dot_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_assemble_ram_issue/NC_001318.1_gene_661/realign.ann.tested.ft.count.ft.tred.dot";
    
    Assembler assembler;
    Graph gp;
    assembler.read_ann_results(ann_file);
    read_dot_file(gp, dot_file);
    assembler.assemble(gp, dot_file + ".assembled");
}

TEST_CASE("test multithread ann", "[hide]")
{
    string encode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_ann_multithread/realign.encode.rdim";
    string recode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_ann_multithread/realign.recode";
    string m5_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_ann_multithread/realign.m5";
    string out_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_ann_multithread/realign.ann";
    
    Assembler assembler;
    assembler.ann_clust_recode_multithread(recode_file, recode_file + ".ref", encode_file, m5_file,
                                           10, 0.2, 0.8, 25, 50, 0.5, false, true, 1, false, 4);
    
    vector<int64_t> idx;
    assembler.print_rl_ann_clust(out_file + ".raw", true);
    assembler.find_nccontigs(idx);
    assembler.print_rl_ann_clust(out_file + ".igda_tmp", true, idx);
    string cmd = "sort -u -s -k2,2n -k3,3n -k1,1 " + out_file + ".igda_tmp" + " > " + out_file;
    cout << cmd << endl; system(cmd.c_str());
}

TEST_CASE("test test_correct_contigs_pairwise_multiple", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_correct_contigs_pairwise/results/realign.ann.tested.ft.count";
    string recode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_correct_contigs_pairwise/results/realign.recode";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_correct_contigs_pairwise/results/Borrelia_burgdorferi_N40_chr.fna";

    Assembler assembler;
    assembler.load_homo_blocks(ref_file);
    assembler.set_null_betadist(1.332824, 89.04769);
    assembler.test_contigs_pairwise(ann_file, recode_file, ann_file + ".ft", 20, 3, 10, 20);
   
}

TEST_CASE("test test_contigs (new algorithm)", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_correct_contigs_pairwise/results_test_contigs/realign.ann.count";
    string recode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_correct_contigs_pairwise/results_test_contigs/realign.recode";
    string align_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_correct_contigs_pairwise/results_test_contigs/realign_readrange.m5";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_correct_contigs_pairwise/results_test_contigs/Borrelia_burgdorferi_N40_chr.fna";
    double alpha = 1.332824;
    double beta = 89.04769;
    double min_logbf = 20;
    int max_loci = 10;
    double min_rr = 20;
    
    Assembler assembler;
    
    cout << "load ann" << endl;
    assembler.read_ann_results(ann_file);
    
    cout << "load recode_data" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    cout << "load recode_ref_data" << endl;
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_file + ".ref");
    
    cout << "load reads_range" << endl;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    cout << "load ref_file" << endl;
    assembler.load_homo_blocks(ref_file);
    
    cout << "test contigs" << endl;
    assembler.set_null_betadist(alpha, beta);
    assembler.test_contigs(recode_data, recode_ref_data, reads_range);
    assembler.print_rl_ann_clust(ann_file + ".tested", true);
    
    cout << "filter contigs" << endl;
    assembler.filter_ann(ann_file + ".tested", min_logbf, max_loci, min_rr);
    
}

TEST_CASE("test find_nccontigs (new overlap rule)", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_new_overlap_rule/realign.ann.tested.ft.count.ft";
    string ann_nc_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_new_overlap_rule/realign.ann.tested.ft.count.ft.nc.j07";
    
    Assembler assembler;
    assembler.read_ann_results(ann_file);
    vector<int64_t> idx;
    assembler.find_nccontigs(idx, 0.5, 0.7);
    
    assembler.print_rl_ann_clust(ann_nc_file + ".igda_tmp", true, idx);
    string cmd = "sort -u -s -k2,2n -k3,3n -k1,1 " + ann_nc_file + ".igda_tmp" + " > " + ann_nc_file;
    cout << cmd << endl; system(cmd.c_str());
    cmd = "rm -f " + ann_nc_file+ ".igda_tmp";
    cout << cmd << endl; system(cmd.c_str());
}

TEST_CASE("test tred (new overlap rule)", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_new_overlap_rule/realign.ann.tested.ft.count.ft.nc.j07";
    
    Assembler assembler;
    Graph gp_bgl;
    //assembler.ann_to_graph(gp_bgl, ann_file, 0.5, 0.5, 2, false);
    assembler.ann_to_graph(gp_bgl, ann_file, 0.5, 0.5, 2, true);
    
    ofstream fs_graph(ann_file + ".dot");
    boost::write_graphviz(fs_graph, gp_bgl);
    fs_graph.close();
    
    // transitive reduction
    IGDA_Graph gp;
    load_igda_graph_from_file(gp, ann_file + ".dot", ann_file);
    
    IGDA_Graph gp_tred;
    igda_tred(gp, gp_tred);
    
    save_igda_graph_to_file(gp_tred, ann_file + ".tred.dot");
    
}

TEST_CASE("test Assemble:cut_overhanged_contigs", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_cut_overhanged_contigs/realign.ann.tested.ft.count.ft.assembled.count.nc.ft";
    string ann_cut_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_cut_overhanged_contigs/realign.ann.tested.ft.count.ft.assembled.count.nc.ft.cut";
    Assembler assembler;
    assembler.cut_overhanged_contigs(ann_file, ann_cut_file);
    //assembler.read_ann_results(ann_file);
}

TEST_CASE("test tred or that minimal length or minimal number of SNVs intead of and", "[hide]")
{
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred_or/unassembled.ann";
    Assembler assembler;
    
    // overlap
    Graph gp_bgl;
    assembler.ann_to_graph(gp_bgl, ann_file, 0.5, 0.5, 0.7, true);
    
    ofstream fs_graph(ann_file + ".dot");
    boost::write_graphviz(fs_graph, gp_bgl);
    fs_graph.close();
    
    // transitive reduction
    IGDA_Graph gp;
    load_igda_graph_from_file(gp, ann_file + ".dot", ann_file);
    
    IGDA_Graph gp_tred;
    igda_tred(gp, gp_tred);
    
    save_igda_graph_to_file(gp_tred, ann_file + ".tred.dot");
    
    
}

