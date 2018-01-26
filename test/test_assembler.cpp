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
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string out_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.jaccard.5000";
    Assembler assembler;
    assembler.jaccard_index(encode_file, m5_file, out_file);
}

TEST_CASE("test assembler::assemble_core()")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    
    vector<vector<int> > encode_data;
    vector<ReadRange> reads_range;
    loadencodedata(encode_data, encode_file);
    loadreadsrange(reads_range, m5_file);
    
    vector<vector<int> > centroid;
    vector<ReadRange> centroid_range;
    vector<int> n_idx_on;
    vector<vector<int> > idx_on;
    Assembler assembler;
    assembler.assemble_core(encode_data, reads_range, centroid, centroid_range, idx_on, n_idx_on, 20, 200, 100);
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.centroid");
    for (int i=0; i<(int)centroid.size(); ++i){
        fs_outfile << centroid[i] << endl;
    }
    fs_outfile.close();
    
    open_outfile(fs_outfile, "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.n_idx_on");
    fs_outfile << n_idx_on << endl;
    fs_outfile.close();
}

TEST_CASE("test assembler::check_contained_reads()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    
    vector<vector<int> > encode_data;
    vector<ReadRange> reads_range;
    loadencodedata(encode_data, encode_file);
    loadreadsrange(reads_range, m5_file);
    
    Assembler assembler;
    vector<bool> is_contained = assembler.check_contained_reads(encode_data, reads_range, 200);
    for (int i=0; i<(int)is_contained.size(); ++i)
        if (!is_contained[i])
            cout << i <<',';
    cout << endl;
}


TEST_CASE("test assembler::mat_fac_rank_1()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    
    vector<vector<int> > encode_data;
    vector<ReadRange> reads_range;
    loadencodedata(encode_data, encode_file);
    loadreadsrange(reads_range, m5_file);
    
    vector<int> idx_on; vector<int> idx_off;
    Assembler assembler;
    ReadRange centroid_range = reads_range[23];
    vector<int> centroid = encode_data[23];
    cout << "centroid_range : [" << centroid_range.first<< "," << centroid_range.second << "]" << endl;
    cout << "centroid : " << centroid << endl;
    int n_iter = assembler.mat_fac_rank_1(encode_data, reads_range, centroid, centroid_range,
                                          idx_on, idx_off, 20, 200, 100);
    cout << "new centroid : " << centroid << endl;
    cout << "idx_on : " << idx_on << endl;
    cout << "number of iteration : " << n_iter << endl;
}



TEST_CASE("test assembler::mat_fac_rank_1_core()", "[hide]")
{
    string encode_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string m5_file = "../results/dforest/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    
    vector<vector<int> > encode_data;
    vector<ReadRange> reads_range;
    loadencodedata(encode_data, encode_file);
    loadreadsrange(reads_range, m5_file);
    
    vector<int> idx_on; vector<int> idx_off;
    Assembler assembler;
    ReadRange centroid_range = reads_range[23];
    vector<int> centroid = encode_data[23];
    cout << "centroid_range : [" << centroid_range.first<< "," << centroid_range.second << "]" << endl;
    cout << "centroid : " << centroid << endl;
    assembler.mat_fac_rank_1_core(encode_data, reads_range, centroid, centroid_range, idx_on, idx_off, 200);
    cout << "new centroid : " << centroid << endl;
    cout << "idx_on : " << idx_on << endl;
}








