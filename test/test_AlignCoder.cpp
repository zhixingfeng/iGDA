//
//  test_AlignCoder.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/modules.h"

TEST_CASE("test AlignCoderSNV")
{
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    p_aligncoder->encode("../data/MSSA_61_forward.m5", "../results/MSSA_61_forward_encode_snv.txt");
}
