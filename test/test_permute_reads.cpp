//
//  test_permute_reads.cpp
//  iGDA
//
//  Created by Zhixing Feng on 1/28/20.
//  Copyright © 2020 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../src/misc/permute_reads.h"
TEST_CASE("test gen_binom", "[hide]")
{
    pcg32 rng(1363);
    //vector<int> cand = {0,1,2,3};
    //vector<double> prob = {0.003428571, 0.014040816, 0.006312925, 0.000000000};
    
    vector<int> sub_cand = {0,1};
    vector<double> sub_prob = {0.9762177, 0.02378231};
    
    
    for (int i = 0; i < 10000; ++i){
        int is_sub = gen_binom(sub_prob, sub_cand, rng);
        if (is_sub == 1){
            vector<int> cand = {0,1,2,3};
            vector<double> prob = {0.003428571, 0.014040816, 0.006312925, 0.000000000};
            int r_num = gen_binom(prob, cand, rng);
            cout << r_num << ",";
        }
    }
    
    
}

TEST_CASE("test permute_encodefile", "[hide]")
{
    string m5_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.m5";
    string pu_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.pileup";
    string out_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.encode.permuted";
    permute_encodefile(m5_file, pu_file, out_file);
}

TEST_CASE("test get_condprob_threshold", "[hide]")
{
    string dforest_permuted_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.dforest.permuted";
    string pu_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.pileup";
    get_condprob_threshold(dforest_permuted_file, pu_file);
}
