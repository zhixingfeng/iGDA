//
//  test_misc.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/10/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//
#include "../include/catch.hpp"
#include "../src/misc/misc.h"


vector<CondFreq> getcondfreq(vector<int> &candidates, vector<vector<int> > &encode_data,
                             vector<ReadRange> &reads_range, vector<int> &pu,
                             vector<double> &p0, string mode="snv")
{
    if (candidates.back() > pu.size())
        throw runtime_error("candidates.back() > pu.size()");
    if (candidates.back() > p0.size())
        throw runtime_error("candidates.back() > p0.size()");

    // get marginal frequency for each candidate
    vector<int> freq_mar(candidates.size(), -1);
    for (int i=0; i<(int)candidates.size(); i++)
        freq_mar[i] = pu[candidates[i]-1];
    
    // get p0 for candidates
    vector<double> cand_p0(candidates.size(), -1);
    for (int i=0; i<(int)candidates.size(); i++)
        cand_p0[i] = p0[candidates[i]-1];
    
    // get binary in the subspace spaned by candidates
    vector<vector<int> > submat= getsubspace(candidates, encode_data, (int)pu.size(), reads_range);
    
    // get conditional freq
    vector<CondFreq> condfreq(candidates.size(),CondFreq());
    for (int j=0; j<(int)submat[0].size(); j++){
        // get index of 1 and -1(missing value) at locus j
        vector<int> idx_1;
        vector<int> idx_1n;
        for (int i=0; i<(int)submat.size(); i++){
            if (submat[i][j]==1)
                idx_1.push_back(i);
            if (submat[i][j]==-1)
                idx_1n.push_back(i);
        }
        if (idx_1.size()==0) continue;
        
        // get conditional frquency with locus k
        for (int k=0; k<(int)submat[0].size(); k++){
            if (k==j) continue;
            CondFreq cur_condfreq(0,freq_mar[k]);
            for (int i=0; i<(int) idx_1.size(); i++)
                if (submat[idx_1[i]][k]==1)
                    cur_condfreq.x++;
            for (int i=0; i<(int) idx_1n.size(); i++)
                if (submat[idx_1n[i]][k]==1)
                    cur_condfreq.n--;
            if (cur_condfreq.x > cur_condfreq.n) throw runtime_error("x>n");
            if (cur_condfreq.n > 0){
                cur_condfreq.p = cur_condfreq.x / cur_condfreq.n;
                if (cur_condfreq.p > p0[j]){
                    cur_condfreq.log_lr = binom_log_lr(cur_condfreq.x, cur_condfreq.n, p0[j]);
                }else{
                    cur_condfreq.log_lr = 0;
                }
            }
            if (cur_condfreq.log_lr > condfreq[j].log_lr){
                condfreq[j] = cur_condfreq;
                condfreq[j].idx = k;
            }
        }
    }
    
    // douoble check marginal frquency
    /*for (int j=0; j<(int)submat[0].size(); j++){
        int cur_mar = 0;
        for (int i=0; i<(int)submat.size(); i++)
            if (submat[i][j]==1)
                cur_mar++;
        if (cur_mar != freq_mar[j])
            throw runtime_error("cur_mar != freq_mar[j]");
    }*/
    return condfreq;
}

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

TEST_CASE("Test cmpreads","[hide]"){
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    string out_file = "../results/B_10_cons_encode_snv_cmpreads.txt";
    
    //vector<vector<int> > encode_data;
    //loadencodedata(encode_data, encode_file);
    
    //vector<ReadRange> reads_range;
    //loadreadsrange(reads_range, align_file);
    
    cmpreads(encode_file, align_file, out_file);
    
}

TEST_CASE("Test pileup", "[hide]")
{
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
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

TEST_CASE("Test condfrq, getcondfreq")
{
    cout <<"Test condfrq, getcondfreq" << endl;
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    
    int candidates_[] = {143,557,629,811,848,1073,1297,1455,1491,1729,1742,1790,2091,2118,2576,3000,4000};
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
        fs_outfile_condfreq << condfreq[i].log_lr << '\t' << condfreq[i].idx + 1 << endl;
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




