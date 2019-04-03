//
//  test_dforeststxxl.cpp
//  iGDA
//
//  Created by Zhixing Feng on 2018/5/26.
//  Copyright © 2018年 Zhixing Feng. All rights reserved.
//


#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/dforest/dforestsnvmax.h"
#include "../src/modules/dforest/dforestsnvstxxl.h"

/*TEST_CASE("test dforest", "[hide]")
{
    string align_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string encode_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.5000";
    string cmpreads_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.cmpreads.5000";
    string out_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.dforest.5000";

    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVMax forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results/dforeststxxl/tmp" , 8, 5, 1);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}

TEST_CASE("test dforeststxxl", "[hide]")
{
    string align_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string encode_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.5000";
    string cmpreads_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.cmpreads.5000";
    string out_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.dforest.stxxl.5000";
    
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVSTXXL forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results/dforeststxxl/tmp" , 8, 5, 1, 0, true);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}


TEST_CASE("test dforeststxxl (whole data)", "[hide]")
{
    string align_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5";
    string encode_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode";
    string cmpreads_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.cmpreads";
    string out_file = "../results/dforeststxxl/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.dforest.stxxl";
    
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVSTXXL forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results/dforeststxxl/tmp" , 8, 5, 1, 0, true);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}


TEST_CASE("test dforeststxxl (use sam file)", "[hide]")
{
    string align_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.sam";
    string encode_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.sam.encode";
    string cmpreads_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.sam.cmpreads";
    string out_file = "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.5000.sam.dforest";
    
    AlignReaderSam alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVSTXXL forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results/dforeststxxl/tmp" , 8, 5, 1, 0, true);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}*/


TEST_CASE("test dforeststxxl report progresss", "[hide]")
{
    string align_file = "../data/loman/even/realign.m5";
    string encode_file = "../data/loman/even/realign.encode";
    string cmpreads_file = "../data/loman/even/realign.cmpreads";
    string ref_file = "../data/loman/Listeria_monocytogenes_complete_genome_masked.fasta";
    string out_file = "../data/loman/igda_even/realign.dforest";
    string tmp_dir = "../data/loman/igda_even/tmp";
    
    AlignReaderSam alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVSTXXL forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    int start_time= (int)clock();
    
    ptr_forest->load_homo_blocks(ref_file);
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, tmp_dir, 15, 10000, 1, 0, 0.75, 25, false);
   
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}

