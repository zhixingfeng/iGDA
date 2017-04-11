//
//  test_errormodel.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/10.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//
#include "../include/catch.hpp"
#include "../src/modules/errormodel/errormodelsnv.h"

TEST_CASE("test ErrorModelSNV::learn()")
{
    string align_file = "../data/B_10_cons_forward.m5";
    ErrorModelSNV model;
    model.learn(align_file, "../results/B_10_cons");
}
