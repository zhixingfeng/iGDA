//
//  test_sclust.cpp
//  iGDA
//
//  Created by Zhixing Feng on 7/17/17.
//  Copyright (c) 2017 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/sclust/sclust.h"

#include <ctime>

TEST_CASE("test sclust::run()")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads_topn.bin";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n1.sclust";
    
    SClust sclust;
    sclust.run(encode_file, align_file, cmpreads_file, out_file, "./", 15, 0, 0, 0, 1);
}
