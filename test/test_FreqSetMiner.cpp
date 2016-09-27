//
//  test_FreqSetMiner.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/9/26.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//


#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/modules.h"


TEST_CASE("Test FreqSetMiner"){
    string cmpreadsfile = "../results/B_10_cons_encode_snv_cmpreads.txt";
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    
    FreqSetMinerSNV FreqSetMinerSNV_obj;
    FreqSetMiner *ptr_FreqSetMiner = &FreqSetMinerSNV_obj;
    
    ptr_FreqSetMiner->mapEncodetoCmpReads(encode_file, cmpreadsfile);
}
