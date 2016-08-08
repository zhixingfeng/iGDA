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

TEST_CASE("test AlignCoderSNV::encode")
{
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    p_aligncoder->encode("../data/MSSA_61_forward.m5", "../results/MSSA_61_forward_encode_snv.txt");
}

TEST_CASE("test AlignCoderSNV::decode")
{
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;

    REQUIRE(p_aligncoder->decode(2).first==1); REQUIRE(p_aligncoder->decode(2).second=='C');
    REQUIRE(p_aligncoder->decode(3).first==1); REQUIRE(p_aligncoder->decode(3).second=='G');
    REQUIRE(p_aligncoder->decode(4).first==1); REQUIRE(p_aligncoder->decode(4).second=='T');
    REQUIRE(p_aligncoder->decode(5).first==2); REQUIRE(p_aligncoder->decode(5).second=='A');
    REQUIRE(p_aligncoder->decode(6).first==2); REQUIRE(p_aligncoder->decode(6).second=='C');
    REQUIRE(p_aligncoder->decode(7).first==2); REQUIRE(p_aligncoder->decode(7).second=='G');
    REQUIRE(p_aligncoder->decode(8).first==2); REQUIRE(p_aligncoder->decode(8).second=='T');
    
    REQUIRE(p_aligncoder->decode(3894).first==974); REQUIRE(p_aligncoder->decode(3894).second=='C');
    REQUIRE(p_aligncoder->decode(372).first==93);REQUIRE(p_aligncoder->decode(372).second=='T');
    
}
