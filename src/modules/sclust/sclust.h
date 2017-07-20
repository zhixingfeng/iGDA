//
//  sclust.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/7/13.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef sclust_h
#define sclust_h

#include "../../../include/headers.h"
#include "../../misc/misc.h"
#include <thread>
#include <mutex>



class SClust
{
public:
    SClust(){}
    virtual ~SClust(){}
    
    void run(string encode_file, string align_file, string cmpreads_file, 
             string out_file, string tmp_dir, int max_cand_size, int min_ratio, 
             int min_count, int min_cvg, int n_thread);
    
protected:
    bool run_thread(string cmpreads_file, string out_file, int max_cand_size, 
                    int min_ratio, int min_count, int min_cvg);
    
    void count_freq(unordered_set<uint32_t> &pattern, int32_t &nreads_cover_all,
                    const vector<int> &cand_loci, vector<int32_t> &temp_id_var, 
                    vector<int32_t> &temp_id_read, vector<int32_t> &temp_count_var);
    
    void test_pattern(unordered_set<uint32_t> &pattern, int32_t nreads_cover_all, vector<int32_t> &temp_count_var,
                      int min_ratio, int min_count, vector<uint32_t> &rl_pattern, vector<double> &rl_ratio, vector<int> &rl_count);
    
    void print_pattern(FILE *p_outfile, const vector<int> &cand_loci, vector<uint32_t> &rl_pattern,
                       vector<double> &rl_ratio, vector<int> &rl_count, int32_t nreads_cover_all);
    
    // legacy
    void print_freq(FILE *p_outfile, const vector<int> &cand_loci, unordered_set<uint32_t> &pattern,
                    int32_t nreads_cover_all, vector<int32_t> &temp_count_var);
    
protected:
    // pileup variants and reads
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_read;
    int64_t nreads;
    
    // bit shift vector
    vector<uint32_t> bit_shift;
};





#endif /* sclust_h */
