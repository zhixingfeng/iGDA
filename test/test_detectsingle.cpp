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


TEST_CASE("test DetectSingleSNV::loadcontext", "[hide]")
{
    DetectSingleSNV detectsinglesnv;
    DetectSingle *p_detectsingle = &detectsinglesnv;
    p_detectsingle->loadcontexteffect("../results/encode_from_sam/NCTC3000.context");
    p_detectsingle->savecontexteffect("../results/encode_from_sam/NCTC3000.context.reconstructed_mincvg_500");
}

TEST_CASE("test DetectSingleSNV::detect", "[hide]")
{
    string context_file = "../results/encode_from_sam/NCTC3000.context";
    string pileup_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.m5.pileup";
    string out_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.m5.pileup.detectsingle";
    DetectSingleSNV detectsinglesnv;
    DetectSingle *p_detectsingle = &detectsinglesnv;
    p_detectsingle->loadcontexteffect(context_file);
    p_detectsingle->detect(pileup_file, out_file);
    
}


TEST_CASE("test DetectSingleSNV::detect (pt)", "[hide]")
{
    string context_file = "../results/pt/NCTC3000.context";
    string pileup_file = "../results/pt/align_to_consensus.pileup";
    string out_file = "../results/pt/align_to_consensus.pileup.detectsingle";
    DetectSingleSNV detectsinglesnv;
    DetectSingle *p_detectsingle = &detectsinglesnv;
    p_detectsingle->loadcontexteffect(context_file);
    p_detectsingle->detect(pileup_file, out_file, log(0.0000001), 0.00001);
    
}

TEST_CASE("test DetectSingleSNV::detect (HBV)")
{
    string context_file = "../results/LoFreq_HBV/NCTC3000.context";
    string pileup_file = "../results/LoFreq_HBV/igda_result_large/realign.pileup";
    string out_file = "../results/LoFreq_HBV/igda_result_large/realign.detectsingle";
    DetectSingleSNV detectsinglesnv;
    DetectSingle *p_detectsingle = &detectsinglesnv;
    p_detectsingle->loadcontexteffect(context_file);
    p_detectsingle->detect(pileup_file, out_file, log(0.0000001), 0.00001);
    
}



