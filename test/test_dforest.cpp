//
//  test_dforest.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/12/2.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/alignreader/alignreaderm5.h"
#include "../src/modules/aligncoder/aligncodersnv.h"
#include "../src/modules/dforest/dforest.h"

TEST_CASE("test DForest::build_tree()")
{
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    string align_file = "../data/B_10_cons.m5";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForest forest(&alignreader, &aligncoder);
    forest.call_pileup_var(encode_file);
    forest.call_pileup_reads(align_file);
    
    
    vector<int> cand_loci({143,557,629,703,819});
    vector<Result> rl;
    vector<int> temp_vec_var(forest.get_n_reads(), -1);
    vector<int> temp_vec_read(forest.get_n_reads(), -1);

    forest.build_tree(cand_loci, rl, temp_vec_var, temp_vec_read, 10, 10);
    
}

