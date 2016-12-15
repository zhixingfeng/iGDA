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
    
    bool run(string encode_file, string align_file, string cmpreads_file, string a_out_file, int min_reads, int max_depth);
    
    void build_tree(const vector<int> &cand_loci, long long &counter, vector<long long> &temp_vec_var, vector<long long> &temp_vec_read, int min_reads, int max_depth);
};


#endif /* dforestsnv_h */
