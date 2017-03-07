//
//  dforestsnv.hpp
//  iGDA
//
//  Created by Zhixing Feng on 16/12/6.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef dforestsnv_h
#define dforestsnv_h

#include "./dforest.h"

class DForestSNV : public DForest
{
public:
    
    DForestSNV() : DForest() {}
    DForestSNV(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder) : DForest(a_p_alignreader, a_p_aligncoder) {}
    
    virtual ~DForestSNV(){}
    
    bool run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread=1);
    
    void build_tree(FILE * p_cmpreads_file, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth);
    
protected:
    bool run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth);
};


#endif /* dforestsnv_h */
