//
//  dforestsnvstxxl.hpp
//  iGDA
//
//  Created by Zhixing Feng on 10/31/20.
//  Copyright © 2020 Zhixing Feng. All rights reserved.
//

#ifndef dforestsnvstxxl_h
#define dforestsnvstxxl_h

#include "./dforest.h"
#include <thread>
#include <mutex>
#include <stxxl/vector>

typedef stxxl::VECTOR_GENERATOR<int, 4, 8, 8*4>::result stxxl_vector_int;

// define stxxl based vector
class stxxl_v_int
{
public:
    stxxl_v_int() : shift(0), dat(NULL), vec_size(0) {}
    stxxl_v_int(const stxxl_vector_int * _dat) : shift(0), dat(_dat), vec_size(0) {}
    virtual ~stxxl_v_int(){}
    
public:
    const int &operator [] (size_t j) {
        if (this->dat == NULL) throw runtime_error("stxxl_v_int: dat is NULL");
        if (j >= this->vec_size) throw runtime_error("j >= this->vec_size");
        return (*dat)[shift + j];
    }
    size_t & size(){return vec_size;}
public:
    size_t shift;
    //const stxxl::vector<int> *dat;
    const stxxl_vector_int *dat;
    size_t vec_size;
};

// define stxxl based vector of vector
class stxxl_vv_int
{
public:
    stxxl_vv_int(){}
    virtual ~stxxl_vv_int(){}
    
public:
    const stxxl_v_int &operator [] (size_t i) { return this->dat_vec[i]; }
    void pileup_encode(string encode_file);
    void pileup_m5(string m5_file);
    void print_dat_vec(string outfile);
private:
    vector<stxxl_v_int> dat_vec;
    //stxxl::vector<int> dat;
    stxxl_vector_int dat;
};


// define DForestSNVSTXXL
class DForestSNVSTXXL : public DForest
{
public:
    
    DForestSNVSTXXL() : DForest() {}
    DForestSNVSTXXL(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder) : DForest(a_p_alignreader, a_p_aligncoder) {}
    
    virtual ~DForestSNVSTXXL(){}
    
    bool run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread=1, double minfreq=0, double maxfreq=1, int min_homo_block_dist = 15, bool isinter=false);
    
    void build_tree(ofstream &fs_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read,
                    vector<double> &p_y_x_archive, unordered_set<int64_t> &idx_mod, int focal_locus, int min_reads, int max_depth, double minfreq, bool isinter=false);
    
    void save_result(string out_file, double minfreq);
    void save_result_all(string out_file, double minfreq);
    
protected:
    
    bool run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq, double maxfreq, int min_homo_block_dist = 15, bool isinter=false);
    void build_index(stxxl_vector_type &cmpreads_index, const stxxl_vector_type_int &cmpreads_data);
    
protected:
    vector<DforestResult> result;
    vector<DforestResult> result_all;
    
    int min_homo_block_dist;
    
    std::mutex thread_locker;
};
#endif /* dforestsnvstxxl_h */
