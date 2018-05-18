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

TEST_CASE("Test seqopt, getrevcomp"){
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

TEST_CASE("Test cmpreads_topn (read data from memory or stxxl)", "[hide]"){
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    string out_txtfile = "../results/B_10_cons_cmpreads_topn.norange.txt.stxxl";
    
    // load encode_data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_data
    AlignReaderM5 AlignReaderM5_obj;
    stxxl::vector<Align> align_data;
    AlignReaderM5_obj.read(align_file, align_data);

    // run cmpreads_topn
    stxxl::vector<vector<int> > cmpreads_data;
    clock_t t_begin = clock();
    cmpreads_topn(encode_data, align_data, cmpreads_data, 10, 0, true, false, false);
    clock_t t_end = clock();

    // print results
    ofstream fs_outfile; open_outfile(fs_outfile, out_txtfile);
    for (int i=0; i<(int)cmpreads_data.size(); ++i)
        fs_outfile << cmpreads_data[i] << "," << endl;
    fs_outfile.close();
    
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
    stxxl::vector<ReadMatch> cmpreadsdiff_data;
    loadcmpreadsdiff(cmpreadsdiff_data, cmpreads_diff_file);
    
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000.reconstruct.txt");
    for (int i = 0; i < (int)cmpreadsdiff_data.size(); ++i){
        fs_outfile << cmpreadsdiff_data[i].read_id << '\t' << cmpreadsdiff_data[i].start << '\t' << cmpreadsdiff_data[i].end << '\t';
        fs_outfile << cmpreadsdiff_data[i].matches << '\t' << cmpreadsdiff_data[i].diff << endl;
    }
    fs_outfile.close();
}

TEST_CASE("test readmatch_compare")
{
    string cmpreads_diff_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000";
    stxxl::vector<ReadMatch> cmpreadsdiff_data;
    loadcmpreadsdiff(cmpreadsdiff_data, cmpreads_diff_file);
    
    stxxl::vector<vector<ReadMatch> > cmpreadsdiff_data_grouped;
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

/*TEST_CASE("test build_index_cmpreadsdiff")
{
    string cmpreads_diff_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000";
    stxxl::vector<ReadMatch> cmpreadsdiff_data;
    loadcmpreadsdiff(cmpreadsdiff_data, cmpreads_diff_file);
    
    stxxl::vector<vector<size_t> > cmpreadsdiff_data_index;
    build_index_cmpreadsdiff(cmpreadsdiff_data, cmpreadsdiff_data_index);
    
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000.reconstruct.fromindex.txt");
    for (size_t i = 0; i < cmpreadsdiff_data_index.size(); ++i){
        for (size_t j = 0; j < cmpreadsdiff_data_index[i].size(); ++j){
            size_t idx = cmpreadsdiff_data_index[i][j];
            fs_outfile << cmpreadsdiff_data[idx].read_id << '\t' << cmpreadsdiff_data[idx].start << '\t' << cmpreadsdiff_data[idx].end << '\t';
            fs_outfile << cmpreadsdiff_data[idx].matches << '\t' << cmpreadsdiff_data[idx].diff << endl;
        }
    }
    
}*/







