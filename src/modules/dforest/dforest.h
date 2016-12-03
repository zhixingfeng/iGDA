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
#include "../../misc/pileup.h"
#include "../alignreader/alignreaderm5.h"
#include "../aligncoder/aligncodersnv.h"

// define final result.
struct Result{
    
    double bf;
    double prop;
    vector<int> link_loci;
    
};

class DForest {
    
public:
    DForest(){p_alignreader = NULL; p_aligncoder = NULL;}
    DForest(AlignReader *a_p_alignreader, AlignCoder *a_p_aligncoder){
        p_alignreader = a_p_alignreader;
        p_aligncoder = a_p_aligncoder;
        p_aligncoder->setAlignReader(p_alignreader);
    }
    virtual ~DForest(){}
    
    inline void call_pileup_var(string encode_file){
        pu_var = pileup_var(encode_file, n_reads);
    }
    inline void call_pileup_reads(string align_file, char format = 'm'){
        long int cur_n_reads;
        pu_read = pileup_reads(align_file, cur_n_reads, format);
        if (cur_n_reads != n_reads)
            throw runtime_error("number of reads in align_file and encode_file are different");
    }
    inline long int get_n_reads(){return n_reads;}
    
    inline void setAlignReader(AlignReader *a_p_alignreader){p_alignreader = a_p_alignreader;}
    inline void setAlignCoder(AlignCoder *a_p_aligncoder){
        p_aligncoder = a_p_aligncoder;
        p_aligncoder->setAlignReader(p_alignreader);
    }
    
    inline bool run(string cmpreads_file, string out_file, int min_reads = 10, int max_depth=3){
        vector<int> temp_vec_var(n_reads, -1);
        vector<int> temp_vec_read(n_reads, -1);
        return false;
    }
    
protected:
    AlignReader *p_alignreader;
    AlignCoder *p_aligncoder;
    
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_read;
    
    long int n_reads;
    
    string out_file;
    
public:
    inline void build_tree(const vector<int> &cand_loci, int locus_y, Result &rl, int min_reads, int max_depth, vector<int> &temp_vec_var, vector<int> & temp_vec_read)
    {
        // calculate joint frequency of variants
        
        // calculate marginal frequency of predictors
        
    }
    
};





#endif /* dforest_h */
