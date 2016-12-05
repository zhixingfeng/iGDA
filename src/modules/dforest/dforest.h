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
    inline void build_tree(const vector<int> &cand_loci, vector<Result> &rl, vector<int> &temp_vec_var, vector<int> & temp_vec_read, int min_reads, int max_depth)
    {
        // each of the locus in cand_loci is used as response y
        for (int i = 0; i < cand_loci.size(); i++){
            // NOTE: cand_loci[i] is response y. Let fill in temp_vec_var and temp_vec_read
            // by reponse y.
            int y_locus = cand_loci[i];
            int y_read_locus = int (y_locus / 4);
            
            for (int j = 0; j < pu_var[y_locus].size(); j++)
                temp_vec_var[pu_var[y_locus][j]] = y_locus;
            
            for (int j = 0; j < pu_read[y_read_locus].size(); j++)
                temp_vec_read[pu_read[y_read_locus][j]] = y_read_locus;
            
            // calculate joint frequency of variants
            vector<int> p_y_x(cand_loci.size(), -1);
            for (int j = 0; j < cand_loci.size(); j++){
                // avoid self comparison. p_y_x will be reused, so give it -1 instead of skipping
                // if we can not get a meaning value;
                if (j == i){
                    p_y_x[j] = -1;
                    continue;
                }
                
                // calculate p_y_x by filling temp_vec_var
                for (int k = 0; k < pu_var[cand_loci[j]].size(); j++){
                    //if ()
                }
                    
            }
            
        }
        
        
    }
    
};





#endif /* dforest_h */
