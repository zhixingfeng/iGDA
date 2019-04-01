//
//  dforestsnvstxxl.hpp
//  iGDA
//
//  Created by Zhixing Feng on 2018/5/26.
//  Copyright © 2018年 Zhixing Feng. All rights reserved.
//

#ifndef dforestsnvstxxl_hpp
#define dforestsnvstxxl_hpp

#include "./dforest.h"



class DForestSNVSTXXL : public DForest
{
public:
    
    DForestSNVSTXXL() : DForest() {}
    DForestSNVSTXXL(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder) : DForest(a_p_alignreader, a_p_aligncoder) {}
    
    virtual ~DForestSNVSTXXL(){}
    
    bool run(const vector<vector<int> > &encode_data, const stxxl::vector<Align> &align_data,
             const stxxl::vector<vector<int> > &cmpreads_data, int min_reads, int max_depth,
             int n_thread=1, double minfreq=0);
    
    bool run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread=1, double minfreq=0, double maxfreq=1, bool isinter=false);
    
    void build_tree(ofstream &fs_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq, bool isinter=false);
    
    // get_result() not actually used in DForestSNVSTXXL, retained for historical reasons
    inline unordered_map<int, DforestResult> get_result()
    {
        return unordered_map<int, DforestResult>();
    }
    
    void save_result(string out_file, double minfreq);
    void save_result_all(string out_file, double minfreq);
    
protected:
    
    bool run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq, double maxfreq, bool isinter=false);
    void build_index(stxxl_vector_type &cmpreads_index, const stxxl_vector_type_int &cmpreads_data);
    
protected:
    vector<DforestResult> result;
    stxxl::vector<DforestResult> result_all;
        
    int focal_locus;
    vector<double> p_y_x_archive;
    unordered_set<int64_t> idx_mod; // record index of modified p_y_x_archive
};




#endif /* dforestsnvstxxl_hpp */
