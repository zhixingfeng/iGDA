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
#include "../src/modules/dforest/dforestsnvmax.h"
#include "../src/modules/dforest/dforestsnvstl.h"
#include <ctime>
#include <stxxl/vector>
/*TEST_CASE("test DForest::run() (input from memory/stxxl)", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads_topn.txt";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n1.txt.stxxl";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVMax forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    cout << "load encode_data" << endl;
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);

    cout << "load alignment" << endl;
    AlignReaderM5 AlignReaderM5_obj;
    vector<Align> align_data;
    AlignReaderM5_obj.read(align_file, align_data);

    cout << "cmpreads" << endl;
    stxxl_vector_type_int cmpreads_data;
    loadcmpreads(cmpreads_data, cmpreads_file);

    cout << "dforest" << endl;
    int start_time= (int)clock();
    ptr_forest->run(encode_data, align_data, cmpreads_data, 8, 5, 1);
    int stop_time= (int)clock();
    
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
    
    // write results (unordered) to outfile
    ofstream fs_outfile;  open_outfile(fs_outfile, out_file);
    unordered_map<int, DforestResult> result = ptr_forest->get_result();
    for (auto it = result.begin(); it!=result.end(); ++it){
        if (it->second.link_loci.size() > 0 && it->second.p_y_xp >= 0){
            fs_outfile << it->second.focal_locus << '\t' << it->second.bf << '\t'
                        << it->second.p_y_xp << '\t' << it->second.n_y_xp << '\t'
                        << it->second.n_xp << '\t' << it->second.link_loci.size() << '\t';
            for (int j = 0; j < it->second.link_loci.size(); j++)
                fs_outfile << it->second.link_loci[j] << ',';
            fs_outfile << endl;
        }
    }
    fs_outfile.close();
    
}*/


/*TEST_CASE("test DForest::run()", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads_topn.bin";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n1.txt";
    //string out_file = "../results/dummy_cmpreads_out.txt";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    //DForestSNV forestsnv(&alignreader, &aligncoder);
    //DForestSNVFast forestsnv(&alignreader, &aligncoder);
    DForestSNVMax forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results" , 8, 5, 1);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}

TEST_CASE("test DForest::filter()", "[hide]")
{
    string dforest_file = "../results/B_10_cons.dforest";
    string out_file = "../results/B_10_cons.dforest.f0.5.filter";
    
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVMax forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    ptr_forest->filter(dforest_file, out_file, 0.5);
}*/


/*TEST_CASE("test DForest::build_tree()", "[hide]")
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
        ptr_forest->build_tree(cand_loci, temp_vec_var, temp_vec_var_lock, temp_vec_read, temp_vec_read_lock, 8, 10);
    clock_t t_end = clock();
    cout << "time for build_tree : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}
*/


/*TEST_CASE("test DForest::run() (hbv)", "[hide]")
{
    string encode_file = "../results/LoFreq_HBV/igda_result_large/realign.encode";
    string align_file = "../results/LoFreq_HBV/igda_result_large/realign.m5";
    string cmpreads_file = "../results/LoFreq_HBV/igda_result_large/realign.cmpreads";
    string ref_file = "../results/LoFreq_HBV/cat_wild_large.fasta";
    string out_file = "../results/LoFreq_HBV/igda_result_large/realign.dforest.debug.v2";
    
    
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVSTL forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    ptr_forest->load_homo_blocks(ref_file);
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results/LoFreq_HBV/igda_result_large/tmp" , 12, 10000, 1);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}*/



