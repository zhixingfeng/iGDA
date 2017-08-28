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


typedef vector<uint32_t> AssemblyNode;
typedef vector<vector<AssemblyNode> > AssemblyGraph;


class SClust
{
public:
    SClust(){}
    virtual ~SClust(){}
    
    void run(string encode_file, string align_file, string cmpreads_file, 
             string out_file, string tmp_dir, int max_cand_size, double min_logLR, 
             int min_count, int min_cvg, int n_thread);
    
    void eval_pattern(string pattern_file, string true_snp_file, string out_file);
    
    // merge pattern of the same read
    void summary(string sclust_file, string out_file, double min_logLR, int min_count, int min_cvg);
    
    // assembly corrected reads
    void assemble(string summary_file, string out_file, int min_count, int min_cvg);

    
    // split long subspaces
    void split_subspace(string cmpreads_file, string out_file, int max_cand_size);
    
    // merge pattern of the same read (legacy)
    void summary_old(string sclust_file, string out_file, double min_ratio, double min_logLR, 
                 int min_count, int min_cvg);

    
protected:
    bool run_thread(string cmpreads_file, string out_file, int max_cand_size, 
                    int min_logLR, int min_count, int min_cvg);
    
    // count frequency of patterns
    void count_freq(unordered_set<uint32_t> &pattern, int32_t &nreads_cover_all,
                    const vector<int> &cand_loci, vector<int32_t> &temp_id_var, 
                    vector<int32_t> &temp_id_read, vector<int32_t> &temp_count_var);
    
    // pairwise test indepence 
    void pairwise_test(vector<vector<double> > &pairwise_logLR, const vector<int> &cand_loci,
                       vector<int32_t> &temp_id_var, vector<int32_t> &temp_id_read);
    
    void test_pattern(const unordered_set<uint32_t> &pattern, int32_t nreads_cover_all, const vector<int> &cand_loci,
                      const vector<int32_t> &temp_count_var, int min_count, double min_logLR,
                      vector<double> &rl_logLR, vector<int> &rl_count);
    void print_pattern(FILE *p_outfile, const int read_id, const int start, const int end, const vector<int> &cand_loci,
                       vector<double> &rl_logLR, vector<int> &rl_count, int32_t nreads_cover_all);

    
    /*void test_pattern(const vector<vector<double> > &pairwise_logLR, const unordered_set<uint32_t> &pattern,
                      int32_t nreads_cover_all, const vector<int32_t> &temp_count_var, 
                      int min_count, double min_logLR,
                      vector<uint32_t> &rl_pattern, vector<double> &rl_logLR, vector<int> &rl_count);*/

    //void print_pattern(FILE *p_outfile, const int read_id, const vector<int> &cand_loci, vector<uint32_t> &rl_pattern,
    //                       vector<double> &rl_logLR, vector<int> &rl_count, int32_t nreads_cover_all);
    
    // test significance of patterns (legacy)
    /*void test_pattern_old(unordered_set<uint32_t> &pattern, int32_t nreads_cover_all, vector<int32_t> &temp_count_var,
                      int min_ratio, int min_count, vector<uint32_t> &rl_pattern, 
                      vector<double> &rl_logLR, vector<double> &rl_ratio, vector<int> &rl_count);
    
    // output pattern (legacy)
    void print_pattern_old(FILE *p_outfile, const int read_id, const vector<int> &cand_loci, vector<uint32_t> &rl_pattern,
                       vector<double> &rl_logLR, vector<double> &rl_ratio, vector<int> &rl_count, 
                       int32_t nreads_cover_all);*/
    
    
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
    
    inline int32_t get_max_count(unordered_set<uint32_t> &pattern, vector<int32_t> &temp_count_var)
    {
        int32_t max_count = 0;
        for (auto it=pattern.begin(); it!=pattern.end(); ++it){
            if (temp_count_var[*it] > max_count && bitcount(*it)>=2)
                max_count = temp_count_var[*it];
        }
        return max_count;
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
