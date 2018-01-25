//
//  assemble.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef assemble_h
#define assemble_h

#include "../../../include/headers.h"
#include "../../misc/io.h"
#include "../../misc/basic.h"
#include "../aligncoder/aligncodersnv.h"

class Assembler 
{
public:
    Assembler(){}
    virtual ~Assembler(){}
    
public:
    void get_variants(string dforest_file, string out_file, double min_condprob);
    void reduce_dim(string encode_file, string var_file, string out_file);
    void dist(string encode_file, string align_file, string out_file);
    void dist_rdim(string encode_file, string align_file, string var_file, string out_file);
    void jaccard_index(string encode_file, string align_file, string out_file);
    
    
    void assemble_core(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                       vector<vector<int> > &centroid, vector<ReadRange> &centroid_range, vector<int> &n_idx_on,
                       int min_idx_on, int min_overlap, int max_iter);
    
    // check if a read is contained by another one
    vector<bool> check_contained_reads(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                                       int min_overlap, bool rm_empty_centroid = true);
    
    // rank 1 matrix facterization, return number of iteration.
    int mat_fac_rank_1(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                        vector<int> &centroid, const ReadRange &centroid_range,
                        vector<int> &idx_on, vector<int> &idx_off, int min_idx_on, int min_overlap, int max_iter);
    
    void mat_fac_rank_1_core(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                             vector<int> &centroid, const ReadRange &centroid_range, 
                             vector<int> &idx_on, vector<int> &idx_off, int min_overlap);
}
;

#endif /* assemble_h */


