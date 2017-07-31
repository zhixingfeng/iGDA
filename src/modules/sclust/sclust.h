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
    
    void eval_pattern(string pattern_file, string true_snp_file, string out_file);
    
    // merge similar patterns
    void summary(string pattern_file, string out_file, int min_overlap);

protected:
    bool run_thread(string cmpreads_file, string out_file, int max_cand_size, 
                    int min_ratio, int min_count, int min_cvg);
    
    // count frequency of patterns
    void count_freq(unordered_set<uint32_t> &pattern, int32_t &nreads_cover_all,
                    const vector<int> &cand_loci, vector<int32_t> &temp_id_var, 
                    vector<int32_t> &temp_id_read, vector<int32_t> &temp_count_var);
    
    // test significance of patterns
    void test_pattern(unordered_set<uint32_t> &pattern, int32_t nreads_cover_all, vector<int32_t> &temp_count_var,
                      int min_ratio, int min_count, vector<uint32_t> &rl_pattern, 
                      vector<double> &rl_logLR, vector<double> &rl_ratio, vector<int> &rl_count);
    
    // output pattern
    void print_pattern(FILE *p_outfile, const vector<int> &cand_loci, vector<uint32_t> &rl_pattern,
                       vector<double> &rl_logLR, vector<double> &rl_ratio, vector<int> &rl_count, 
                       int32_t nreads_cover_all);
    
    
//protected:
public:
    inline double cal_logLR(double n_11, double n_10, double n_01, double N)
    {
        double n_00 = N - n_11 - n_10 - n_01;
        if ( N*n_11 > (n_11 + n_10)*(n_11 + n_01) ){
            return cal_logL_H1(n_11, n_10, n_01, n_00, N) - cal_logL_H0(n_11, n_10, n_01, n_00, N);
        }else{
            return 0;
        }
    }
    inline double cal_logL_H0(double n_11, double n_10, double n_01, double n_00, double N)
    {
        double n_1x = n_11 + n_10; double n_0x = n_01 + n_00;
        double n_x1 = n_11 + n_01; double n_x0 = n_10 + n_00;
        return cal_nlogn(n_1x) + cal_nlogn(n_0x) + cal_nlogn(n_x1) + cal_nlogn(n_x0) - 2*cal_nlogn(N);
    }
    inline double cal_logL_H1(double n_11, double n_10, double n_01, double n_00, double N)
    {
        return cal_nlogn(n_11) + cal_nlogn(n_10) + cal_nlogn(n_01) + cal_nlogn(n_00) - cal_nlogn(N);
    }
    inline double cal_nlogn(double x)
    {
        if (int(x) == 0)
            return 0;
        return x*log(x);
    }
    
protected:
    // pileup variants and reads
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_read;
    int64_t nreads;
    
    // bit shift vector
    vector<uint32_t> bit_shift;
};





#endif /* sclust_h */
