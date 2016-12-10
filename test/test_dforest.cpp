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
#include "../src/modules/dforest/dforestsnv.h"


TEST_CASE("test DForest::run()")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads.bin";
    string out_file = "../results/B_10_cons_out.txt";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNV forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    ptr_forest->run(align_file, encode_file, cmpreads_file, out_file, 8, 5);
}
TEST_CASE("test DForest::build_tree()", "[hide]")
{
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNV forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    ptr_forest->call_pileup_var(encode_file);
    ptr_forest->call_pileup_reads(align_file);
    
    vector<int> cand_loci({142, 556, 628, 702, 818, 1001, 1003, 1050, 10000});
    vector<vector<Result> > rl(ptr_forest->get_pileup_var().size(), vector<Result>());
    vector<int> temp_vec_var(ptr_forest->get_n_reads(), -1);
    vector<int> temp_vec_var_lock(ptr_forest->get_n_reads(), -1);
    vector<int> temp_vec_read(ptr_forest->get_n_reads(), -1);
    vector<int> temp_vec_read_lock(ptr_forest->get_n_reads(), -1);

    clock_t t_begin = clock();
    for (int i=0; i<1000; i++)
        ptr_forest->build_tree(cand_loci, rl, temp_vec_var, temp_vec_var_lock, temp_vec_read, temp_vec_read_lock, 8, 10);
    clock_t t_end = clock();
    cout << "time for build_tree : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}






