//
//  dforest.h
//  iGDA
//
//  Created by Zhixing Feng on 16/12/1.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//
//  To implement data-based forest method for bump hunting, we demond large amount of computation as we have to train p*n^2 models,
//  where p is average shared variants of two reads and n is number of reads. We need fast yet low memory cost ways to train models.
//  For large datasets, read ID pileup file might not be able to fit into memory, so we use two different way to implement: first is
//  the naive way, and the second is prefix compression (memory efficient but might be slower).

//  Input is pileup data, i.e Pileup variant and pileup read ID.
//  Parameters: minimal reads supporting a partition, maximal depth of a tree.


#ifndef dforest_h
#define dforest_h


#include "../../../include/headers.h"
#include "../../misc/misc.h"
#include "../alignreader/alignreaderm5.h"
#include "../aligncoder/aligncodersnv.h"
#include <thread>
#include <mutex>
#include <stxxl.h>

// define final result. (locus here means coded locus)
struct DforestResult{
    DforestResult():focal_locus(-1),bf(0),p_y_xp(0),n_y_xp(0),n_xp(0){}
    int focal_locus; 
    double bf;
    double p_y_xp;
    int n_y_xp;
    int n_xp;
    vector<int> link_loci;
    
};

// define DForest model. Note that all the locus stored in the files are 1-based
// pileup read ID is 0-based

class DForest {
    
public:
    DForest(){p_alignreader = NULL; p_aligncoder = NULL;}
    DForest(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder){
        p_alignreader = a_p_alignreader;
        p_aligncoder = a_p_aligncoder;
        p_aligncoder->setAlignReader(p_alignreader);
    }
    virtual ~DForest(){}
    
    inline void call_pileup_var(const vector<vector<int> > &encode_data){
        pu_var = pileup_var(encode_data, n_reads);
    }
    inline void call_pileup_reads(const stxxl::vector<Align> &align_data, char format = 'm'){
        int64_t cur_n_reads;
        pu_read = pileup_reads(align_data, cur_n_reads, format);
        if (cur_n_reads != n_reads)
            throw runtime_error("number of reads in align_file and encode_file are different");
    }

    inline void call_pileup_var(string encode_file){
        pu_var = pileup_var(encode_file, n_reads);
    }
    inline void call_pileup_reads(string align_file, char format = 'm'){
        int64_t cur_n_reads;
        pu_read = pileup_reads(align_file, cur_n_reads, format);
        if (cur_n_reads != n_reads)
            throw runtime_error("number of reads in align_file and encode_file are different");
    }
    
    inline int64_t get_n_reads(){return n_reads;}
    inline vector<vector<int> > get_pileup_var(){return pu_var;}
    inline vector<vector<int> > get_pileup_reads(){return pu_read;}
    
    inline void setAlignReader(AlignReader *a_p_alignreader){p_alignreader = a_p_alignreader;}
    inline void setAlignCoder(AlignCoder *a_p_aligncoder){
        p_aligncoder = a_p_aligncoder;
        p_aligncoder->setAlignReader(p_alignreader);
    }
    
    virtual bool run(const vector<vector<int> > &encode_data, const stxxl::vector<Align> &align_data,
                     const stxxl::vector<vector<int> > &cmpreads_data, int min_reads, int max_depth,
                     int n_thread=1, double minfreq=0)=0;
    
    virtual bool run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread=1, double minfreq=0)=0;
        
    virtual void build_tree(FILE * p_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq) = 0;
    
    virtual unordered_map<int, DforestResult> get_result()=0;
    void filter(string dforest_file, string out_file, double minfreq);
    
protected:
    AlignReader *p_alignreader;
    AlignCoder *p_aligncoder;
    
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_read;
    
    int64_t n_reads;
    
    
    
    
    
};


#endif /* dforest_h */
