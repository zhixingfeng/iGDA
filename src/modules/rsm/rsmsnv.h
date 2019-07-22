//
//  rsmsnv.hpp
//  iGDA
//
//  Created by Zhixing Feng on 7/22/19.
//  Copyright © 2019 Zhixing Feng. All rights reserved.
//

#ifndef rsmsnv_h
#define rsmsnv_h


#include "../dforest/dforest.h"

class RSMsnv{
public:
    RSMsnv(){}
    RSMsnv(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder){
        p_alignreader = a_p_alignreader;
        p_aligncoder = a_p_aligncoder;
    }
    virtual ~RSMsnv(){}
    
    inline void load_homo_blocks(const string fasta_file)
    {
        this->homo_blocks = get_homo_blocks(fasta_file);
    }
    
    bool run(string encode_file, string encode_ref_file, string align_file, string cmpreads_file, string out_file, int min_reads, int max_depth, int n_thread=1, double minfreq=0, double maxfreq=1, int min_homo_block_dist = 15, bool isinter=false);
    void build_tree(ofstream &fs_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq, bool isinter=false);
    void save_result(string out_file, double minfreq);
    
protected:
    bool run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq, double maxfreq, int min_homo_block_dist = 15, bool isinter=false);
    void build_index(stxxl_vector_type &cmpreads_index, const stxxl_vector_type_int &cmpreads_data);
    
protected:
    AlignReader *p_alignreader;
    AlignCoder *p_aligncoder;
    
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_var_ref;
    vector<vector<int> > pu_read;
    
    int64_t n_reads;
    
    vector<int64_t> homo_blocks;
    
    vector<DforestResult> result;
    vector<DforestResult> result_all;
    
    int focal_locus;
    vector<double> p_y_x_archive;
    unordered_set<int64_t> idx_mod; // record index of modified p_y_x_archive
    
    int min_homo_block_dist;
    
    
};


#endif /* rsmsnv_hpp */
