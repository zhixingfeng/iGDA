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
    pcg32 rng(17363);
    vector<int> cand = {0,1,2,3};
    vector<double> prob = {0.003428571, 0.014040816, 0.006312925, 0.000000000};
    
    for (int i = 0; i < 10000; ++i)
        cout << gen_binom(prob, cand, rng) << ",";
    
    
}

TEST_CASE("test permute_encodefile")
{
    string encode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.encode";
    string pu_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.pileup";
    string out_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_permutation/detect/clpA_1/realign.encode.permuted";
    permute_encodefile(encode_file, pu_file, out_file);
}

