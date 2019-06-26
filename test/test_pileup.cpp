//
//  test_pileup.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/11/30.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/misc.h"

TEST_CASE("test pileup_qv()", "[hide]")
{
    string sam_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/pileup/SAMEA4916110_split/NZ_UWXV01000001.1_forward.sam";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/pileup/GCF_900608245.1_kpneu039_genomic_NZ_UWXV01000001.1.fna";
    
    vector<vector<pair<int64_t, double> > > pu_qv = pileup_qv(sam_file, ref_file);
    
    int x = 1;
    
}

TEST_CASE("test pileup_var() (input from memory/stxxl)", "[hide]")
{
    clock_t t_begin = clock();
    string encode_file = "../results/B_10_cons.encode";
    int64_t n_reads;
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    vector<vector<int> > pu_var = pileup_var(encode_data, n_reads);
    clock_t t_end = clock();
    cout << "time for pileup_var() : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    ofstream fs_out; open_outfile(fs_out, "../results/B_10_cons.pileup_var.txt.mem.stxxl");
    for (int i = 0; i < pu_var.size(); i++){
        for (int j = 0; j < pu_var[i].size(); j++){
            fs_out << pu_var[i][j] << '\t';
        }
        fs_out << endl;
    }
    fs_out.close();
}


TEST_CASE("test pileup_var() (input from file)", "[hide]")
{
    clock_t t_begin = clock();
    string encode_file = "../results/B_10_cons.encode";
    int64_t n_reads;
    vector<vector<int> > pu_var = pileup_var(encode_file, n_reads);
    clock_t t_end = clock();
    cout << "time for pileup_var() : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    ofstream fs_out; open_outfile(fs_out, "../results/B_10_cons.pileup_var.txt");
    for (int i = 0; i < pu_var.size(); i++){
        for (int j = 0; j < pu_var[i].size(); j++){
            fs_out << pu_var[i][j] << '\t';
        }
        fs_out << endl;
    }
    fs_out.close();
}

TEST_CASE("test pileup_reads() (input from memory/stxxl)", "[hide]")
{
    clock_t t_begin = clock();
    string align_file = "../data/B_10_cons.m5";
    AlignReaderM5 AlignReaderM5_obj;
    vector<Align> align_data;
    AlignReaderM5_obj.read(align_file, align_data);
    
    int64_t n_reads;
    vector<vector<int> > pu_read = pileup_reads(align_data, n_reads, 'm');
    clock_t t_end = clock();
    cout << "time for pileup_reads() : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    ofstream fs_out; open_outfile(fs_out, "../results/B_10_cons_pileup_reads.txt.mem.stxxl");
    for (int i = 0; i < pu_read.size(); i++){
        for (int j = 0; j < pu_read[i].size(); j++){
            fs_out << pu_read[i][j] << '\t';
        }
        fs_out << endl;
    }
    
    fs_out.close();
}

TEST_CASE("test pileup_reads() (input from file)", "[hide]")
{
    clock_t t_begin = clock();
    string align_file = "../data/B_10_cons.m5";
    int64_t n_reads;
    vector<vector<int> > pu_read = pileup_reads(align_file, n_reads, 'm');
    clock_t t_end = clock();
    cout << "time for pileup_reads() : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    ofstream fs_out; open_outfile(fs_out, "../results/B_10_cons_pileup_reads.txt");
    for (int i = 0; i < pu_read.size(); i++){
        for (int j = 0; j < pu_read[i].size(); j++){
            fs_out << pu_read[i][j] << '\t';
        }
        fs_out << endl;
    }
    
    fs_out.close();
}

TEST_CASE("test print_pileup()", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    int64_t n_reads;
    vector<vector<int> > pu_var = pileup_var(encode_file, n_reads);
    vector<vector<int> > pu_reads = pileup_reads(align_file, n_reads);
    
    print_pileup(pu_var, "../results/B_10_cons.var.pilup");
    print_pileup(pu_reads, "../results/B_10_cons.reads.pilup");
}


TEST_CASE("test pileup_reads() removing deletions", "[hide]")
{
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.1000";
    
    int64_t n_reads;
    vector<vector<int> > pu_reads = pileup_reads(align_file, n_reads, false);
    print_pileup(pu_reads, "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.reads.pileup.1000");
    
    vector<vector<int> > pu_reads_rm_del = pileup_reads(align_file, n_reads, true);
    print_pileup(pu_reads_rm_del, "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.reads.pileup.rm_del.1000");
}

TEST_CASE("test filter_pileup_var()", "[hide]")
{
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.1000";
    string recode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.1000";
    int64_t n_reads;
    int64_t n_reads_var;
    vector<vector<int> > pu_reads_rm_del = pileup_reads(align_file, n_reads, true);
    vector<vector<int> > pu_var = pileup_var(recode_file, n_reads_var);
    if (n_reads_var != n_reads)
        throw runtime_error("n_reads_var != n_reads");
    
    vector<vector<int> > pu_var_ft = filter_pileup_var(pu_var, pu_reads_rm_del, n_reads);
    print_pileup(pu_var_ft, "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.var_ft.pileup.1000");
}







