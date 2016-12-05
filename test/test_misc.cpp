//
//  test_misc.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/10/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//
#include "../include/catch.hpp"
#include "../src/misc/misc.h"



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

TEST_CASE("Test cmpreads"){
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    string out_file = "../results/B_10_cons_encode_snv_cmpreads_array_method_rm_single.txt";
    
    clock_t t_begin = clock();
    cmpreads(encode_file, align_file, out_file, 0);
    clock_t t_end = clock();
    cout << "time for compare reads : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
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

TEST_CASE("Test lgamma")
{
    REQUIRE(lgamma(1)==0);REQUIRE(lgamma(2)==0);REQUIRE(lgamma(4)==Approx(1.791759e+00).epsilon(0.00001));
    REQUIRE(lgamma(64)==Approx(201.0093).epsilon(0.00001));REQUIRE(lgamma(1024)==Approx(6071.28).epsilon(0.00001));
    REQUIRE(lgamma(65536)==Approx(661276.9).epsilon(0.00001));REQUIRE(lgamma(1048576)==Approx(13487768).epsilon(0.00001));
}






