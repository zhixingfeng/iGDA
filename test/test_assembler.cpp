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


TEST_CASE("test assembler::ann_clust")
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
    assembler.ann_clust_recode_legacy(recode_file, recode_ref_file, m5_file, var_file, 12, 0.2, 0.8, 20, 50, 0.5);
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
    assembler.ann_clust_recode_legacy(recode_file, recode_ref_file, m5_file, var_file, 12, 0.2, 0.8, 20, 50, 0.5);

}

TEST_CASE("test assembler::assign_reads_to_contigs()", "[hide]")
{
    //string ann_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann";
    //string recode_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode";
    //string m5_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.m5";
    string ann_file = "../results/pt_rpoBC/igda/realign.ann.tested.ft";
    string recode_file = "../results/pt_rpoBC/igda/realign.recode";
    string m5_file = "../results/pt_rpoBC/igda/realign.m5";
    
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


TEST_CASE("test assembler::ann_to_graph()", "[hide]")
{
    Assembler assembler;
    Graph gp, gp_red;
    assembler.ann_to_graph(gp, "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft");
    
    // transitive reduction using BGL, but it seems incorrect
    igda_transitive_reduction(gp, gp_red);
    print_graph(gp);
    cout << "transitive reduction" << endl;
    print_graph(gp_red);
    
    // transitive reduction using graphviz
    ofstream fs_temp_graph("../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft.dot");
    boost::write_graphviz(fs_temp_graph, gp);
    
    string cmd = "igda_tred ../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft.dot ";
    cmd = cmd + "> ../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft.transitive_reduction.dot";
    cout << cmd << endl;
    system(cmd.c_str());
    
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


TEST_CASE("test assembler::test_contigs() (hbv)", "[hide]")
{
    string ref_file = "../results/LoFreq_HBV/cat_wild_large.fasta";
    string ann_file = "../results/LoFreq_HBV/igda_result_large/realign.ann";
    string recode_file = "../results/LoFreq_HBV/igda_result_large/realign.recode";
    string recode_ref_file = "../results/LoFreq_HBV/igda_result_large/realign.recode.ref";
    string m5_file = "../results/LoFreq_HBV/igda_result_large/realign.m5";
    string out_file = "../results/LoFreq_HBV/igda_result_large/realign.ann.tested";
    
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
    string recode_file = "../results/LoFreq_HBV/igda_result_large/3l/realign.recode";
    string recode_ref_file = "../results/LoFreq_HBV/igda_result_large/3l/realign.recode.ref";
    string encode_file = "../results/LoFreq_HBV/igda_result_large/3l/realign.encode.rdim";
    string m5_file = "../results/LoFreq_HBV/igda_result_large/3l/realign.m5";
    string var_file = "../results/LoFreq_HBV/igda_result_large/3l/realign.var";
    string out_file = "../results/LoFreq_HBV/igda_result_large/3l/align_to_consensus_trim.ann.debug";
    
    Assembler assembler;
    
    assembler.ann_clust_recode(recode_file, recode_ref_file, encode_file, m5_file, var_file, 10, 0.2, 0.8, 12, 30, 0.5);
    
    vector<int64_t> idx;
    assembler.find_nccontigs(idx);
    assembler.print_rl_ann_clust(out_file, true, idx);
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





