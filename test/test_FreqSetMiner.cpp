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


TEST_CASE("Test FreqSetMiner::detectVariantsSingle", "[hide]"){
    string cmpreadsfile = "../results/B_10_cons_encode_snv_cmpreads.txt";
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    FreqSetMinerSNV FreqSetMinerSNV_obj;
    FreqSetMiner *ptr_FreqSetMiner = &FreqSetMinerSNV_obj;
    
    AlignCoderSNV aligncodersnv;
    ptr_FreqSetMiner->setAlignCoder(&aligncodersnv);
    vector<double> p_values = ptr_FreqSetMiner->detectVariantsSingle(encode_file, align_file);
    
    ofstream fs_outfile; open_outfile(fs_outfile, "../results/B_10_cons_detect_single_pvalues.txt");
    for (int i=0; i<(int)p_values.size(); i++){
        fs_outfile << p_values[i] << endl;
    }
    fs_outfile.close();
    
}

TEST_CASE("Test FreqSetMiner", "[hide]"){
    string cmpreadsfile = "../results/B_10_cons_encode_snv_cmpreads.txt";
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    FreqSetMinerSNV FreqSetMinerSNV_obj;
    FreqSetMiner *ptr_FreqSetMiner = &FreqSetMinerSNV_obj;
    
    AlignCoderSNV aligncodersnv;
    ptr_FreqSetMiner->setAlignCoder(&aligncodersnv);
    vector<int> variable_loci = ptr_FreqSetMiner->detectVariantsCoarse(encode_file, align_file, cmpreadsfile, 0.01);
    
    
}

