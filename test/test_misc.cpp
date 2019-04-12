//
//  test_misc.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/10/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//
#include "../include/catch.hpp"
#include "../src/misc/misc.h"
#include <thread>
#include <stxxl.h>

TEST_CASE("Test m5tofa", "[hide]"){
    m5tofa("../data/MSSA_61_forward.m5", "../results/MSSA_61_forward.fa");
}

TEST_CASE("Test seqopt, getrevcomp", "[hide]"){
    string originseq = "CTAGTCNN-N--CGTAGT-CGT-CGATGCTGTAGCTA";
    string revseq = getrevcomp(originseq);
    REQUIRE(revseq == "TAGCTACAGCATCG-ACG-ACTACG--N-NNGACTAG");
    originseq = "CGT-CGTNN--CG";
    revseq = getrevcomp(originseq);
    REQUIRE(revseq == "CG--NNACG-ACG");
}

TEST_CASE("Test loadencodedata", "[hide]"){
    string encode_file = "../data/SM_263.code";
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
}


TEST_CASE("Test cmpreads_topn_diff with read_id", "[hide]"){
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string align_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_binfile = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.cmpreads.diff";
    string out_txtfile = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.cmpreads.diff.txt";
    
    clock_t t_begin = clock();
    cmpreads_topn_diff(encode_file, align_file, out_binfile, 10, 0);
    clock_t t_end = clock();
    cout << "time for compare reads (binary output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;

    cmpreads_bin2txt(out_binfile, out_txtfile, false, true);
}


TEST_CASE("Test cmpreads_topn with read_id is_condprob = true","[hide]"){
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    string out_txtfile = "../results/B_10_cons_cmpreads_topn_readid_condprob.txt";
    string out_binfile = "../results/B_10_cons_cmpreads_topn_readid_condprob.bin";
    
    clock_t t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_txtfile, 10, 0, true, false, true, true);
    clock_t t_end = clock();
    cout << "time for compare reads (text output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_binfile, 10, 0, true, true, true, true);
    t_end = clock();
    cout << "time for compare reads (binary output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
}

TEST_CASE("Test cmpreads_topn with read_id","[hide]"){
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    string out_txtfile = "../results/B_10_cons_cmpreads_topn_readid.txt";
    string out_binfile = "../results/B_10_cons_cmpreads_topn_readid.bin";
    
    clock_t t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_txtfile, 10, 0, true, false, true, false);
    clock_t t_end = clock();
    cout << "time for compare reads (text output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_binfile, 10, 0, true, true, true, false);
    t_end = clock();
    cout << "time for compare reads (binary output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    
}

TEST_CASE("Test cmpreads_bin2txt with read_id", "[hide]")
{
    string cmpreads_binfile = "../results/B_10_cons_cmpreads_topn_readid.bin";
    string cmpreads_txtfile = "../results/B_10_cons_cmpreads_topn_readid.bin.totxt";
    cmpreads_bin2txt(cmpreads_binfile, cmpreads_txtfile , true);
    
    //cmpreads_binfile = "../results/B_10_cons_cmpreads_topn_readid_condprob.bin";
    //cmpreads_txtfile = "../results/B_10_cons_cmpreads_topn_readid_condprob.bin.totxt";
    //cmpreads_bin2txt(cmpreads_binfile, cmpreads_txtfile , true);

    
}


TEST_CASE("Test cmpreads_topn", "[hide]"){
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    string out_txtfile = "../results/B_10_cons_cmpreads_topn.norange.txt";
    string out_binfile = "../results/B_10_cons_cmpreads_topn.norange.bin";
    
    clock_t t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_txtfile, 10, 0, true, false, false);
    clock_t t_end = clock();
    cout << "time for compare reads (text output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_binfile, 10, 0, true, true, false);
    t_end = clock();
    cout << "time for compare reads (binary output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}



TEST_CASE("Test cmpreads_bin2txt", "[hide]")
{
    string cmpreads_binfile = "../results/B_10_cons_cmpreads.bin";
    string cmpreads_txtfile = "../results/B_10_cons_cmpreads.bin2txt";
    cmpreads_bin2txt(cmpreads_binfile, cmpreads_txtfile);
}

TEST_CASE("Test cmpreads_txt2bin", "[hide]")
{
    cout << "Test cmpreads_txt2bin";
    string cmpreads_txtfile = "../results/B_10_cons_cmpreads.txt";
    string cmpreads_binfile = "../results/B_10_cons_cmpreads.txt2bin";
    cmpreads_txt2bin(cmpreads_txtfile, cmpreads_binfile);
}

TEST_CASE("Test cmpreads_split", "[hide]")
{
    string cmpreads_binfile = "../results/B_10_cons_cmpreads.bin";
    cmpreads_split(cmpreads_binfile, "../results/B_10_cons_cmpreads.bin.part", 5);
}


TEST_CASE("Test pnorm","[hide]")
{
    cout << pnorm(-10) << endl;
    cout << pnorm(-5) << endl;
    cout << pnorm(-3) << endl;
    cout << pnorm(-1) << endl;
    cout << pnorm(0) << endl;
    cout << pnorm(0.5) << endl;
    cout << pnorm(2.1) << endl;
    cout << pnorm(5) << endl;
    cout << pnorm(10) << endl;
}

TEST_CASE("Test lgamma","[hide]")
{
    REQUIRE(lgamma(1)==0);REQUIRE(lgamma(2)==0);REQUIRE(lgamma(4)==Approx(1.791759e+00).epsilon(0.00001));
    REQUIRE(lgamma(64)==Approx(201.0093).epsilon(0.00001));REQUIRE(lgamma(1024)==Approx(6071.28).epsilon(0.00001));
    REQUIRE(lgamma(65536)==Approx(661276.9).epsilon(0.00001));REQUIRE(lgamma(1048576)==Approx(13487768).epsilon(0.00001));
}

TEST_CASE("Test sort_order", "[hide]")
{
    vector<int> x_int = {7,4,7,3,1,5,6,2,3,6};
    vector<int> idx_int_inc = sort_order(x_int, false);
    vector<int> idx_int_dec = sort_order(x_int, true);
    
    vector<double> x_double = {7.0, 4.0, 7.0, 3.0, 1.0, 5.0, 6.0, 2.5, 3.0, 6.0};
    vector<int> idx_double_inc = sort_order(x_double, false);
    vector<int> idx_double_dec = sort_order(x_double, true);
    
    cout << "x_int: " << x_int << endl;
    cout << "idx_int_inc (increasing): " << idx_int_inc << endl;
    cout << "sorted x_int (increasing): ";
    for (int i=0; i < idx_int_inc.size(); i++)
        cout << x_int[idx_int_inc[i]] << ',';
    cout << endl;
    
    cout << "idx_int_dec (decreasing): " << idx_int_dec << endl;
    cout << "sorted x_int (decreasing): ";
    for (int i=0; i < idx_int_dec.size(); i++)
        cout << x_int[idx_int_dec[i]] << ',';
    cout << endl;
    
    cout << "x_double: " << x_double << endl;
    cout << "idx_double_inc (increasing): " << idx_double_inc << endl;
    cout << "sorted x_double (increasing): ";
    for (int i=0; i < idx_double_inc.size(); i++)
        cout << x_double[idx_double_inc[i]] << ',';
    cout << endl;
    
    cout << "idx_double_dec (decreasing): " << idx_double_dec << endl;
    cout << "sorted x_double (decreasing): ";
    for (int i=0; i < idx_double_dec.size(); i++)
        cout << x_double[idx_double_dec[i]] << ',';
    cout << endl;
    
}

TEST_CASE("Test number of cores", "[hide]")
{
    cout << "cores available: " << thread::hardware_concurrency() << endl;
}

TEST_CASE("Test range of int64_t", "[hide]")
{
    cout << "maximum of int64_t: " << numeric_limits<int64_t>::max() << endl;
}



TEST_CASE("Test bitcount", "[hide]")
{
    cout << bitcount(3) << endl;
}

TEST_CASE("Test slide_win", "[hide]")
{
    vector<int> x = {1,2,3,4,5,6,7,8,9,10};
    int win_size = 2;
    for (int i=0; i<(int)x.size(); ++i){
        vector<int> y;
        slide_win(x, y, i, win_size);
        cout << y << endl;
    }
        
}

TEST_CASE("Test read_fasta", "[hide]")
{
    string ref_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.ref";
    unordered_map<string, string> fasta_data;
    read_fasta(ref_file, fasta_data);
    for (auto it = fasta_data.begin(); it != fasta_data.end(); ++it){
        cout << ">" << it->first << endl;
        cout << it->second << endl;
    }
}


TEST_CASE("Test load_ncread()", "[hide]")
{
    unordered_map<int, int> ncreads = load_ncread("../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected.non_contained.check_follower", 5);
    for (auto it = ncreads.begin(); it != ncreads.end(); ++it)
        cout << it->first << '\t' << it->second << endl;
}


TEST_CASE("test loadcmpreadsdiff", "[hide]")
{
    string cmpreads_diff_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000";
    vector<ReadMatch> cmpreadsdiff_data;
    loadcmpreadsdiff(cmpreadsdiff_data, cmpreads_diff_file);
    
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000.reconstruct.txt");
    for (int i = 0; i < (int)cmpreadsdiff_data.size(); ++i){
        fs_outfile << cmpreadsdiff_data[i].read_id << '\t' << cmpreadsdiff_data[i].start << '\t' << cmpreadsdiff_data[i].end << '\t';
        fs_outfile << cmpreadsdiff_data[i].matches << '\t' << cmpreadsdiff_data[i].diff << endl;
    }
    fs_outfile.close();
}

TEST_CASE("test readmatch_compare", "[hide]")
{
    string cmpreads_diff_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000";
    vector<ReadMatch> cmpreadsdiff_data;
    loadcmpreadsdiff(cmpreadsdiff_data, cmpreads_diff_file);
    
    vector<vector<ReadMatch> > cmpreadsdiff_data_grouped;
    group_cmpreadsdiff(cmpreadsdiff_data, cmpreadsdiff_data_grouped, true);
    
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000.grouped.txt");
    for (size_t i = 0; i < cmpreadsdiff_data_grouped.size(); ++i){
        for (size_t j = 0; j < cmpreadsdiff_data_grouped[i].size(); ++j){
            fs_outfile << cmpreadsdiff_data_grouped[i][j].read_id << '\t' << cmpreadsdiff_data_grouped[i][j].start << '\t' << cmpreadsdiff_data_grouped[i][j].end << '\t';
            fs_outfile << cmpreadsdiff_data_grouped[i][j].matches << '\t' << cmpreadsdiff_data_grouped[i][j].diff << endl;
        }
    }
    fs_outfile.close();
}

TEST_CASE("test get_var_cdf()", "[hide]")
{
    string var_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    vector<int> var_cdf;
    get_var_cdf(var_cdf, var_file, 50000);
    cout << var_cdf << endl;
}

TEST_CASE("test dist_hamming()", "[hide]")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.5000";
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    string var_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    size_t genome_size = get_genome_size(reads_range);
    cout << "genome_size = " << genome_size << endl;
    
    vector<int> var_cdf; get_var_cdf(var_cdf, var_file, genome_size);

    cout << var_cdf << endl;
    
    int i = 0, j = 2;
    vector<bool> temp_array(genome_size*4+3, false);
    double cur_dist = dist_hamming(encode_data[i], encode_data[j], reads_range[i], reads_range[j], var_cdf, temp_array);
    cout << cur_dist << endl;
    
    i = 11; j = 22;
    cur_dist = dist_hamming(encode_data[i], encode_data[j], reads_range[i], reads_range[j], var_cdf, temp_array);
    cout << cur_dist << endl;
    
    i = 4990; j = 4999;
    cur_dist = dist_hamming(encode_data[i], encode_data[j], reads_range[i], reads_range[j], var_cdf, temp_array);
    cout << encode_data[i] << endl;
    cout << encode_data[j] << endl;
    cout << cur_dist << endl;
    
    cout << get_nvar(reads_range[i], reads_range[j], var_cdf) << endl;
    
}


TEST_CASE("test pileup_var (input from encode_data)", "[hide]")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.5000";
    
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // pileup from file
    int64_t n_reads;
    vector<vector<int> > pu_var_file = pileup_var(encode_file, n_reads);
    
    // pileup from data
    vector<vector<int> > pu_var_data = pileup_var(encode_data);
    
    // online pileup
    int pu_size = get_pu_var_size(encode_data);
    vector<vector<int> > pu_var_online(pu_size, vector<int>());
    for (auto i = 0; i < encode_data.size(); ++i)
        pileup_var_online(pu_var_online, encode_data[i], i);
    
    REQUIRE(pu_var_file == pu_var_data);
    REQUIRE(pu_var_file == pu_var_online);
}


TEST_CASE("test pileup_reads_m5 (input from reads_range)", "[hide]")
{
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    // pileup from file
    int64_t n_reads;
    vector<vector<int> > pu_reads_file = pileup_reads_m5(align_file, n_reads, false);
    
    // pileup from data
    vector<vector<int> > pu_reads_data = pileup_reads_m5(reads_range);
    
    // online
    int pu_size = get_pu_read_size(reads_range);
    vector<vector<int> > pu_reads_online(pu_size, vector<int>());
    for (auto i = 0; i < reads_range.size(); ++i)
        pileup_reads_m5_online(pu_reads_online, reads_range[i], i);
    
    REQUIRE(pu_reads_file == pu_reads_data);
    REQUIRE(pu_reads_file == pu_reads_online);
}


TEST_CASE("test pileup_var (input from encode_data) subset", "[hide]")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.5000";
    
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    vector<int> idx = {1,2,5,6,8, 104, 3958};
   
    // pileup from data
    vector<vector<int> > pu_var_data = pileup_var(encode_data, idx);
    
    // online pileup
    int pu_size = get_pu_var_size(encode_data, idx);
    vector<vector<int> > pu_var_online(pu_size, vector<int>());
    for (auto i : idx)
        pileup_var_online(pu_var_online, encode_data[i], i);
    
    REQUIRE(pu_var_data == pu_var_online);
}

TEST_CASE("test pileup_var_count (input from encode_data) subset", "[hide]")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.5000";
    
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    vector<int> idx = {1,2,5,6,8, 104, 3958};
    
    // pileup from data
    vector<vector<int> > pu_var = pileup_var(encode_data, idx);
    vector<int> pu_var_count_ctrl(pu_var.size(),0);
    for (auto i = 0; i<pu_var.size(); ++i)
        pu_var_count_ctrl[i] = (int)pu_var[i].size();
    
    // pileup count
    vector<int> pu_var_count = pileup_var_count(encode_data, idx);
    
    // pileup count online
    int pu_size = get_pu_var_size(encode_data, idx);
    vector<int> pu_var_online_count(pu_size, 0);
    unordered_set<int64_t> mod_idx;
    for (auto i : idx)
        pileup_var_online_count(pu_var_online_count, encode_data[i], mod_idx);

    
    REQUIRE(pu_var_count_ctrl == pu_var_count);
    REQUIRE(pu_var_count_ctrl == pu_var_online_count);
   
}


TEST_CASE("test pileup_reads_m5 (input from reads_range) subset", "[hide]")
{
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    vector<int> idx = {1,2,5,6,8, 104, 3958};
    // pileup from data
    vector<vector<int> > pu_reads_data = pileup_reads_m5(reads_range, idx);
    
    // online
    int pu_size = get_pu_read_size(reads_range, idx);
    vector<vector<int> > pu_reads_online(pu_size, vector<int>());
    for (auto i : idx)
        pileup_reads_m5_online(pu_reads_online, reads_range[i], i);
    
    REQUIRE(pu_reads_data == pu_reads_online);
}

TEST_CASE("test pileup_reads_m5_count (input from reads_range) subset", "[hide]")
{
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    vector<int> idx = {1,2,5,6,8, 104, 3958};
    // pileup from data
    vector<vector<int> > pu_reads = pileup_reads_m5(reads_range, idx);
    vector<int> pu_reads_count_ctrl(pu_reads.size(), 0);
    for (auto i = 0; i < pu_reads.size(); ++i)
        pu_reads_count_ctrl[i] = (int)pu_reads[i].size();
    
    // pileup count
    vector<int> pu_reads_count = pileup_reads_m5_count(reads_range, idx);
    
    // pileup count online
    int pu_size = get_pu_read_size(reads_range, idx);
    vector<int> pu_reads_online_count(pu_size, 0);
    ReadRange mod_range;
    for (auto i : idx)
        pileup_reads_m5_online_count(pu_reads_online_count, reads_range[i], mod_range);

    REQUIRE(pu_reads_count_ctrl == pu_reads_count);
    REQUIRE(pu_reads_count_ctrl == pu_reads_online_count);
}

TEST_CASE("test get_consensus", "[hide]")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.5000";
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    vector<int> idx = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
    vector<int> pu_var_count = pileup_var_count(encode_data, idx);
    vector<int> pu_read_count = pileup_reads_m5_count(reads_range, idx);
    
    ConsensusSeq cons;
    get_consensus(cons, pu_var_count, pu_read_count, 0, 100, 20);
    cout << cons.cons_seq << endl;
}


TEST_CASE("test cmpreads with priority queque", "[hide]")
{
    string align_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string encode_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.5000";
    string cmpreads_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.cmpreads.5000";
    
    cmpreads_topn(encode_file, align_file, cmpreads_file, 20, 0.5);
}


TEST_CASE("test loadreadsrange with sam format", "[hide]")
{
    string m5_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    string sam_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.sam";
    
    cout << "load m5" << endl;
    vector<ReadRange> reads_range_m5;
    loadreadsrange(reads_range_m5, m5_file);
    
    cout << "load sam" << endl;
    vector<ReadRange> reads_range_sam;
    loadreadsrange(reads_range_sam, m5_file);
    
    cout << "get reads_range" << endl;
    string range_m5_file = m5_file + ".range";
    string range_sam_file = sam_file + ".range";
    ofstream fs_out;
    
    open_outfile(fs_out, range_m5_file);
    for (auto cur_range : reads_range_m5)
        fs_out << cur_range.first << ',' << cur_range.second << endl;
    fs_out.close();
    
    open_outfile(fs_out, range_sam_file);
    for (auto cur_range : reads_range_sam)
        fs_out << cur_range.first << ',' << cur_range.second << endl;
    fs_out.close();
    
    string cmd = "sort " + range_m5_file + " | md5";
    system(cmd.c_str());
    
    cmd = "sort " + range_sam_file + " | md5";
    system(cmd.c_str());
}


TEST_CASE("test cmpreads from sam file", "[hide]")
{
    string m5_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.m5";
    string sam_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.sam";
    string encode_m5_file = m5_file + ".encode";
    string encode_sam_file = sam_file + ".encode";
    string cmpreads_m5_file = m5_file + ".cmpreads";
    string cmpreads_sam_file = sam_file + ".cmpreads";
    
    cout << "cmpreads m5 file " << endl;
    cmpreads_topn(encode_m5_file, m5_file, cmpreads_m5_file);
    
    cout << "cmpreads sam file " << endl;
    cmpreads_topn(encode_sam_file, sam_file, cmpreads_sam_file);
    
}

TEST_CASE("test convert sam to m5", "[hide]")
{
    string sam_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.sam";
    string m5_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.m5.fromsam";
    string ref_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.ref.fa";
    
    AlignReaderSam alignreadersam;
    alignreadersam.samtom5(sam_file, ref_file, m5_file);

    string cmd = "cut -f 1,8,9,10,17,18,19,20 -d ' ' ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.m5 | md5";
    system(cmd.c_str());
    
    cmd = "cut -f 1,8,9,10,17,18,19,20 -d ' ' ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.m5.fromsam | md5";
    system(cmd.c_str());
    
}

TEST_CASE("test binom_log_bf", "[hide]")
{
    cout << "binom_log_bf(10, 100, 0.1) = " << binom_log_bf(10, 100, 0.1) << endl;
    cout << "binom_log_bf(1, 100, 0.1) = " << binom_log_bf(1, 100, 0.1) << endl;
    cout << "binom_log_bf(1, 100, 0.5) = " << binom_log_bf(1, 100, 0.5) << endl;
    cout << "binom_log_bf(1, 10, 0.01) = " << binom_log_bf(1, 10, 0.01) << endl;
    cout << "binom_log_bf(1, 20, 0.01) = " << binom_log_bf(1, 20, 0.01) << endl;
    cout << "binom_log_bf(5, 20, 0.01) = " << binom_log_bf(5, 20, 0.01) << endl;
    cout << "binom_log_bf(20, 1000, 0.01) = " << binom_log_bf(20, 1000, 0.01) << endl;
    cout << "binom_log_bf(50, 1000, 0.01) = " << binom_log_bf(50, 1000, 0.01) << endl;
    cout << "binom_log_bf(4, 373, 0.006426349) = " << binom_log_bf(4, 373, 0.006426349) << endl;
}

TEST_CASE("test binom_log_bf with prior of null hypothesis", "[hide]")
{
    cout << "binom_log_bf(1, 20000) = " << binom_log_bf(1, 20000, ALPHA_NULL, BETA_NULL) << endl;
    cout << "binom_log_bf(10, 20000) = " << binom_log_bf(10, 20000, ALPHA_NULL, BETA_NULL) << endl;
    cout << "binom_log_bf(200, 20000) = " << binom_log_bf(200, 20000, ALPHA_NULL, BETA_NULL) << endl;
    cout << "binom_log_bf(1000, 20000) = " << binom_log_bf(1000, 20000, ALPHA_NULL, BETA_NULL) << endl;
    cout << "binom_log_bf(1500, 20000) = " << binom_log_bf(1500, 20000, ALPHA_NULL, BETA_NULL) << endl;
    cout << "binom_log_bf(2000, 20000) = " << binom_log_bf(2000, 20000, ALPHA_NULL, BETA_NULL) << endl;
    cout << "binom_log_bf(3000, 20000) = " << binom_log_bf(3000, 20000, ALPHA_NULL, BETA_NULL) << endl;
    cout << "binom_log_bf(5000, 20000) = " << binom_log_bf(5000, 20000, ALPHA_NULL, BETA_NULL) << endl;
}

TEST_CASE("test cmpreads (fast version)", "[hide]")
{
    string encode_file = "../results/pt/igda_results/align_to_consensus.encode.5000";
    string m5_file = "../results/pt/igda_results/align_to_consensus.m5.5000";
    string cmpreads_file = "../results/pt/igda_results/align_to_consensus.cmpreads.5000";
    
    cmpreads_topn(encode_file, m5_file, cmpreads_file, 20, 0.5);
    
}

TEST_CASE("test sim_jaccard", "[hide]")
{
    string recode_file = "../results/pt_ann_recode_debug/align_to_consensus_trim.recode";
    string m5_file = "../results/pt_ann_recode_debug/align_to_consensus_trim.m5";
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, m5_file);
    
    size_t genome_size = get_genome_size(reads_range);
    
    // create a template to compare reads
    vector<bool> temp_array(genome_size*4+3, false);
    
    int i = 17038;
    int j = 21283;
    double jaccard_index = sim_jaccard(recode_data[i], recode_data[j], reads_range[i], reads_range[j], temp_array);
    
    cout << "recode " << i << ": " << recode_data[i] << endl;
    cout << "range" << i << ": " << reads_range[i].first << ',' << reads_range[i].second << endl;
    
    cout << "recode " << j << ": " << recode_data[j] << endl;
    cout << "range" << j << ": " << reads_range[j].first << ',' << reads_range[j].second << endl;
    
    cout << "jaccard index : " << jaccard_index << endl;

}

TEST_CASE("test get_homo_blocks","[hide]")
{
    
    vector<int64_t> homo_blocks = get_homo_blocks("../results/pt_ann_assign_reads/consensus.fasta");
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/pt_ann_assign_reads/consensus_homo_blocks.txt");
    for (auto i = 0; i < homo_blocks.size(); ++i)
        fs_outfile << homo_blocks[i] << endl;
    fs_outfile.close();
}

TEST_CASE("test prod", "[hide]")
{
    vector<double> x = {1.1, 0.0483, 9.593, 0.1};
    cout << prod(x) << endl;
}

TEST_CASE("test load_varfile", "[hide]")
{
    string var_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.var";
    cout << load_varfile(var_file) << endl;
}

TEST_CASE("cal_locus_specific_mi", "[hide]")
{
    string recode_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode";
    string recode_ref_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ref";
    string var_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.var";
    string m5_file = "../results/pt_ann_assign_reads/align_to_consensus_trim.m5";
    
    vector<vector<int> > recode_data, recode_ref_data;
    loadencodedata(recode_data, recode_file);
    loadencodedata(recode_ref_data, recode_ref_file);
    vector<int> pu_var_count = pileup_var_count(recode_data);
    vector<int> pu_var_ref_count = pileup_var_count(recode_ref_data);
    
    vector<int> var_encode = load_varfile(var_file);
    
    vector<double> weights_11, weights_10, weights_01, weights_00;
    cal_locus_specific_mi(pu_var_count, pu_var_ref_count, var_encode, weights_11, weights_10, weights_01, weights_00);
    
}








