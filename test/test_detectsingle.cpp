//
//  test_detectsingle.cpp
//  iGDA
//
//  Created by Zhixing Feng on 7/5/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../src/modules/detectsingle/detectsingle.h"
#include "../src/modules/detectsingle/detectsinglesnv.h"


TEST_CASE("test loadcontext")
{
    DetectSingleSNV detectsinglesnv;
    detectsinglesnv.loadcontexteffect("../results/encode_from_sam/NCTC3000.context");
    detectsinglesnv.savecontexteffect("../results/encode_from_sam/NCTC3000.context.reconstructed_mincvg_500");
    
}