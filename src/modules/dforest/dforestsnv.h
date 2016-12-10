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
    
    bool run(string align_file, string encode_file, string cmpreads_file, string out_file, int min_reads, int max_depth);
    
    void build_tree(const vector<int> &cand_loci, vector<vector<Result> > &rl, vector<int> &temp_vec_var, vector<int> &temp_vec_var_lock, vector<int> &temp_vec_read, vector<int> &temp_vec_read_lock, int min_reads, int max_depth);
};


#endif /* dforestsnv_h */
