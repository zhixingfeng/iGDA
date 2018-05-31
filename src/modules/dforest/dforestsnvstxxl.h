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
    
    bool run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread=1, double minfreq=0);
    
    void build_tree(FILE * p_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq);
    
    inline unordered_map<int, DforestResult> get_result()
    {
        return result;
    }
    
    void save_result(string out_file, double minfreq);
    
protected:
    
    bool run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq);
    void build_index(vector<stxxl::vector<int64_t> > &cmpreads_index, const stxxl::vector<vector<int> > &cmpreads_data);
    
protected:
    unordered_map<int, DforestResult> result;
    unordered_map<int, DforestResult>::iterator it;
};




#endif /* dforestsnvstxxl_hpp */
