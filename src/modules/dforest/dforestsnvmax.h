//
//  dforestsnvmax.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/6/26.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef dforestsnvmax_h
#define dforestsnvmax_h

#include "./dforest.h"


class DForestSNVMax : public DForest
{
public:
    
    DForestSNVMax() : DForest() {}
    DForestSNVMax(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder) : DForest(a_p_alignreader, a_p_aligncoder) {}
    
    virtual ~DForestSNVMax(){}
    
    bool run(const vector<vector<int> > &encode_data, const stxxl::vector<Align> &align_data,
                     const stxxl::vector<vector<int> > &cmpreads_data, int min_reads, int max_depth,
                     int n_thread=1, double minfreq=0);
    
    bool run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread=1, double minfreq=0, double maxfreq=1, int min_homo_block_dist = 15, bool isinter=false);
    
    void build_tree(ofstream &fs_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq, bool isinter=false);
    
    inline unordered_map<int, DforestResult> get_result()
    {
        return result;
    }
protected:
    bool run_thread_stxxl(const stxxl::vector<vector<int> > &cmpreads_data, int min_reads, int max_depth, double minfreq);
    
    bool run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq, bool isinter=false);

protected:
    unordered_map<int, DforestResult> result;
    unordered_map<int, DforestResult>::iterator it;
};


#endif /* dforestsnvmax_hpp */
