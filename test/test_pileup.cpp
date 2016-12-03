//
//  test_pileup.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/11/30.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/pileup.h"


TEST_CASE("test pileup_var()", "[hide]")
{
    clock_t t_begin = clock();
    string encode_file = "../results/B_10_cons.m5_both_strand_encode_snv.txt";
    long int n_reads;
    vector<vector<int> > pu_var = pileup_var(encode_file, n_reads);
    clock_t t_end = clock();
    cout << "time for pileup_var() : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    ofstream fs_out; open_outfile(fs_out, "../results/B_10_cons.m5_both_strand_encode_snv.pileup_var.txt");
    for (int i = 0; i < pu_var.size(); i++){
        for (int j = 0; j < pu_var[i].size(); j++){
            fs_out << pu_var[i][j] << '\t';
        }
        fs_out << endl;
    }
    fs_out.close();
}

TEST_CASE("test pileup_reads()", "[hide]")
{
    clock_t t_begin = clock();
    string align_file = "../results/B_10_cons_sorted.m5";
    long int n_reads;
    vector<vector<int> > pu_read = pileup_reads(align_file, n_reads, 'm');
    clock_t t_end = clock();
    cout << "time for pileup_reads() : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    ofstream fs_out; open_outfile(fs_out, "../results/B_10_cons_sorted_pileup_reads.txt");
    for (int i = 0; i < pu_read.size(); i++){
        for (int j = 0; j < pu_read[i].size(); j++){
            fs_out << pu_read[i][j] << '\t';
        }
        fs_out << endl;
    }
    
    fs_out.close();
}



