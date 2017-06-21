//
//  test_hclust.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/18.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../src/modules/hclust/hclust.h"

TEST_CASE("test HClust", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string region_file = "../results/B_10_cons.region";
    
    AlignCoderSNV aligncoder;
    HClust hclust(&aligncoder);
    hclust.mask(encode_file, region_file, "../results/B_10_cons_region.clean.encode", false);
    hclust.dist(encode_file, align_file, "../results/B_10_cons.dist", false);
    
}


