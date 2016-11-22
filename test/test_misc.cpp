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

TEST_CASE("Test loadencodedata", "[]"){
    string encode_file = "../data/SM_263.code";
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
}

TEST_CASE("Test cmpreads","[hide]"){
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    string out_file = "../results/B_10_cons_encode_snv_cmpreads.txt";
    
    //vector<ReadRange> reads_range;
    //loadreadsrange(reads_range, align_file);
    
    cmpreads(encode_file, align_file, out_file);
    
}

TEST_CASE("Test pileup", "[hide]")
{
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    //string encode_file = "../data/SM_263.code";
    vector<int> pu = pileup(encode_file);
}

TEST_CASE("Test getcvg", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    vector<int> cvg = getcvg(align_file);
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

TEST_CASE("Test condfrq, getsubspace", "[hide]")
{
    cout <<"Test condfrq, getsubspace" << endl;
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    
    int candidates_[] = {143,557,629,811,848,1073,1297,1455,1491,1729,1742,1790,2091,2118,2576};
    vector<int> candidates (candidates_, candidates_ + sizeof(candidates_)/sizeof(int));
    
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);

    vector<int> pu = pileup(encode_file);
    vector<double> p0(pu.size(), 0.01);
    
    
    vector<vector<int> > submat = getsubspace(candidates, encode_data, (int)pu.size(), reads_range);
    
}

TEST_CASE("Test condfrq, getcondfreq", "[hide]")
{
    cout <<"Test condfrq, getcondfreq" << endl;
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    
    //int candidates_[] = {143,557,629,811,848,1073,1297,1455,1491,1729,1742,1790,2091,2118,2576,3000,4000};
    int candidates_[] = {143,557,703,848,1063,1073,1297,1491,1742,1933,2091,2576,3389,3591,3655,3784,3793,4228,4345,4766,4833,4883,5109,5343,5402,5511,6247,6546,7017,7133,7243,7357,7399,7495,8480,9219,9340,9539,9863,10512,10809,11566,12070,12215,12636,12743,12870,13095,13383,13735,13768,13824,14238,14375,14520,14541,14827,15174,15304,15800,16040,16163,16372,16640,17533,17712,18219,18520,18662,19239,19547,19585,19771,19906,20112,20255,20365,20370,20373,20619,20900,21022,22460,22674,22794,22895,22904,23246};
    vector<int> candidates (candidates_, candidates_ + sizeof(candidates_)/sizeof(int));
    
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    vector<int> pu = pileup(encode_file);
    vector<double> p0(pu.size(), 0.01);
    
    vector<CondFreq> condfreq = getcondfreq(candidates, encode_data, reads_range, pu, p0);
    vector<vector<int> > submat = getsubspace(candidates, encode_data, (int)pu.size(), reads_range);
    
    string outfile_condfreq = "../results/B_10_cons_confreq.txt";
    string outfile_submat = "../results/B_10_cons_submat.txt";
    ofstream fs_outfile_condfreq;
    open_outfile(fs_outfile_condfreq, outfile_condfreq);
    for (int i=0; i<(int)condfreq.size(); i++)
        fs_outfile_condfreq << condfreq[i].log_bf << '\t' << condfreq[i].idx + 1 << endl;
    fs_outfile_condfreq.close();
    
    ofstream fs_outfile_submat;
    open_outfile(fs_outfile_submat, outfile_submat);
    for (int i=0; i<(int)submat.size(); i++){
        for (int j=0; j<(int)submat[i].size(); j++){
            fs_outfile_submat << submat[i][j] << '\t';
        }
        fs_outfile_submat << endl;
    }
    fs_outfile_submat.close();
}




