//
//  test_errormodel.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/10.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//
#include "../include/catch.hpp"
#include "../src/modules/errormodel/errormodelsnv.h"

TEST_CASE("test ErrorModelSNV::learn()", "[hide]")
{
    string align_file = "../data/B_10_cons_forward.m5";
    ErrorModelSNV model;
    model.set_context_size(10, 10);
    model.learn(align_file, "../results/B_10_cons");
}


TEST_CASE("test ErrorModelSNV::merge_all()","[hide]")
{
    ErrorModelSNV model;
    vector<string> context_effect_all_files;
    context_effect_all_files.push_back("../results/contexteffect/ERR718769.context.all");
    context_effect_all_files.push_back("../results/contexteffect/ERR768072.context.all");
    context_effect_all_files.push_back("../results/contexteffect/ERR768076.context.all");
    model.merge_all(context_effect_all_files);
    
}
