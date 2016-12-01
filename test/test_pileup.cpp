//
//  test_pileup.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/11/30.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/pileup.h"


TEST_CASE("test pileup_var()")
{
    clock_t t_begin = clock();
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    vector<vector<int> > pu_var = pileup_var(encode_file);
    clock_t t_end = clock();
    cout << "time for pileup data : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
}

