//
//  assemble.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef assemble_h
#define assemble_h

#include <stxxl.h>
#include "../../../include/headers.h"
#include "../../misc/misc.h"
#include "../aligncoder/aligncodersnv.h"
#include "../src/modules/alignment/alignment.h"
#include "../src/modules/dforest/dforestsnvmax.h"

class Assembler
{
public:
    Assembler():cand_size(5),resampling_size(20),min_count(15),min_condprob(0.15),max_condprob(0.75){}
    virtual ~Assembler(){}
    
public:
    void get_variants(string dforest_file, string out_file, double min_condprob);
    void reduce_dim(string encode_file, string var_file, string out_file);
    void dist(string encode_file, string align_file, string out_file);
    void dist_rdim(string encode_file, string align_file, string var_file, string out_file);
    void jaccard_index(string encode_file, string align_file, string out_file);
    
    inline void call_pileup_var(string encode_file){
        pu_var.clear();
        pu_var = pileup_var(encode_file, n_reads);
    }
    
    inline void call_pileup_reads(string align_file, char format = 'm'){
        pu_read.clear();
        int64_t cur_n_reads;
        pu_read = pileup_reads(align_file, cur_n_reads, format);
        if (cur_n_reads != n_reads)
            throw runtime_error("number of reads in align_file and encode_file are different");
    }
    
    inline void set_par(int cand_size, int resampling_size, int min_count, double min_condprob, double max_condprob)
    {
        cand_size = cand_size;
        resampling_size = resampling_size;
        min_count = min_count;
        min_condprob = min_condprob;
        max_condprob = max_condprob;
    }

    void correct_reads(string encode_file, string align_file, string cmpreads_diff_file, string out_file);
    void run(string encode_file, string align_file, string out_file);
    
    
protected:
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_read;
    int64_t n_reads;
    
    int cand_size;
    int resampling_size;
    int min_count;
    double min_condprob;
    double max_condprob;
    
    vector<vector<int> > adj_mat;
}
;


#endif /* assemble_h */


