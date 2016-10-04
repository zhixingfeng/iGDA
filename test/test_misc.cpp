//
//  test_misc.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/10/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//
#include "../include/catch.hpp"
#include "../src/misc/misc.h"

TEST_CASE("Test m5tofa"){
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

TEST_CASE("Test cmpreads"){
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    //string encode_file = "../data/MSSA_61_forward_encode_snv.txt";
    //string align_file = "../data/MSSA_61_forward.m5";
    string out_file = "../results/B_10_cons_encode_snv_cmpreads.txt";
    
    //vector<vector<int> > encode_data;
    //loadencodedata(encode_data, encode_file);
    
    //vector<ReadRange> reads_range;
    //loadreadsrange(reads_range, align_file);
    
    //cmpreads(encode_file, align_file, out_file);
    
}

TEST_CASE("Test pileup")
{
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    vector<int> pu = pileup(encode_file);
}



