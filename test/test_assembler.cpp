//
//  test_assemble.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/io.h"
#include "../src/misc/basic.h"
#include "../src/modules/assemble/assembler.h"


TEST_CASE("test assembler::get_variants()", "[hide]")
{
    string dforest_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.dforest.uniq.max";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    Assembler assembler;
    assembler.get_variants(dforest_file, out_file, 0.8);
    
}

TEST_CASE("test assembler::reduce_dim()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode";
    string var_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim";
    Assembler assembler;
    assembler.reduce_dim(encode_file, var_file, out_file);
}

TEST_CASE("test assembler::dist()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.dist";
    Assembler assembler;
    assembler.dist(encode_file, m5_file, out_file);
}

TEST_CASE("test assembler::dist_rdim()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string var_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.top20.var";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.dist.withvar";
    Assembler assembler;
    assembler.dist_rdim(encode_file, m5_file, var_file, out_file);
}

TEST_CASE("test assembler::jaccard_index()", "[hide]")
{
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.jaccard.5000";
    Assembler assembler;
    assembler.jaccard_index(encode_file, m5_file, out_file, 0.5);
}

TEST_CASE("test assembler::jaccard_index_min()", "[hide]")
{
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.jaccard_min.5000";
    Assembler assembler;
    assembler.jaccard_index_min(encode_file, m5_file, out_file, 0.1);
}


TEST_CASE("test assembler::correct_reads()", "[hide]")
{
    string align_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string cmpreads_diff_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.cmpreads.diff";
    string out_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected";
    
    Assembler assembler;
    assembler.correct_reads(encode_file, align_file, cmpreads_diff_file, out_file);
}

TEST_CASE("test assembler::check_contained_reads()", "[hide]")
{
    //string align_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    //string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string align_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected.m5";
    string encode_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected";
    vector<vector<int> > encode_data; loadencodedata(encode_data, encode_file);
    vector<ReadRange> reads_range; loadreadsrange(reads_range, align_file);
    
    Assembler assembler;
    vector<int> read_sel_idx = assembler.check_contained_reads(encode_data, reads_range);
    
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected.noncontained_idx");
    fs_outfile << read_sel_idx << endl;
    
    fs_outfile.close();
    cout << read_sel_idx.size() << endl;
    
    select_lines(read_sel_idx, encode_file, encode_file + ".non_contained");
    select_lines(read_sel_idx, align_file, align_file + ".non_contained");
}

TEST_CASE("test assembler::olc()")
{
    string align_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected.m5";
    string encode_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected";
    vector<vector<int> > encode_data; loadencodedata(encode_data, encode_file);
    vector<ReadRange> reads_range; loadreadsrange(reads_range, align_file);
    
    Assembler assembler;
    assembler.olc(encode_file, align_file, "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected.out");
}






