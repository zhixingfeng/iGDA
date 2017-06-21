//
//  dforestsnvfast.h
//  iGDA
//
//  Created by Zhixing Feng on 6/21/17.
//  Copyright (c) 2017 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__dforestsnvfast__
#define __iGDA__dforestsnvfast__

#include "./dforest.h"

class DForestSNVFAST : public DForest
{
public:
    
    DForestSNVFAST() : DForest() {}
    DForestSNVFAST(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder) : DForest(a_p_alignreader, a_p_aligncoder) {}
    
    virtual ~DForestSNVFAST(){}
    
    bool run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread=1, double minfreq=0);
    
    void build_tree(FILE * p_cmpreads_file, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq);
    
protected:
    bool run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq);
    
protected:
    vector<vector<int> > cache_n_y_x;
    vector<vector<int> > cache_n_x;
};

#endif /* defined(__iGDA__dforestsnvfast__) */
