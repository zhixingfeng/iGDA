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

TEST_CASE("Test cmpreads_topn","[hide]"){
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    string out_txtfile = "../results/B_10_cons_cmpreads_topn.txt";
    string out_binfile = "../results/B_10_cons_cmpreads_topn.bin";
    
    clock_t t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_txtfile, 10, 0, true, false);
    clock_t t_end = clock();
    cout << "time for compare reads (text output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    t_begin = clock();
    cmpreads_topn(encode_file, align_file, out_binfile, 10, 0, true, true);
    t_end = clock();
    cout << "time for compare reads (binary output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}


TEST_CASE("Test cmpreads","[hide]"){
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    string out_txtfile = "../results/B_10_cons_cmpreads.txt";
    string out_binfile = "../results/B_10_cons_cmpreads.bin";
    
    clock_t t_begin = clock();
    cmpreads(encode_file, align_file, out_txtfile, 0, true, false);
    clock_t t_end = clock();
    cout << "time for compare reads (text output): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    t_begin = clock();
    cmpreads(encode_file, align_file, out_binfile, 0, true, true);
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

TEST_CASE("Test number of cores")
{
    cout << "cores available: " << thread::hardware_concurrency() << endl;
}

TEST_CASE("Test range of int64_t")
{
    cout << "maximum of int64_t: " << numeric_limits<int64_t>::max() << endl;
}

