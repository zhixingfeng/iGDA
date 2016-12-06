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

// define final result.
struct Result{
    Result():bf(0),p_y_xp(0),n_y_xp(0),n_xp(0){}
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
    inline vector<vector<int> > get_pileup_var(){return pu_var;}
    inline vector<vector<int> > get_pileup_reads(){return pu_read;}
    
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
    inline void build_tree(const vector<int> &cand_loci, vector<vector<Result> > &rl, vector<int> &temp_vec_var, vector<int> &temp_vec_var_lock, vector<int> &temp_vec_read, vector<int> &temp_vec_read_lock, int min_reads, int max_depth)
    {
        // each of the locus in cand_loci is used as response y
        vector<double> p_y_x(cand_loci.size(), -1);
        for (int i = 0; i < cand_loci.size(); i++){
            Result cur_rl;
            // NOTE: cand_loci[i] is response y. Let's fill in temp_vec_var and temp_vec_read
            // by reponse y. pu_var is 4 times larger pu_read because of binary coding so
            // we have to devide 4 to access pu_read
            int y_locus = cand_loci[i];
            int y_read_locus = int (y_locus / 4);
            
            // fill in temp_vec_var and temp_vec_read by response y
            for (int j = 0; j < pu_var[y_locus].size(); j++)
                temp_vec_var[pu_var[y_locus][j]] = y_locus;
            
            for (int j = 0; j < pu_read[y_read_locus].size(); j++)
                temp_vec_read[pu_read[y_read_locus][j]] = y_read_locus;
            
            // calculate joint frequency of neighbor variants
            for (int j = 0; j < cand_loci.size(); j++){
                // avoid self comparison. p_y_x will be reused, so give it -1 instead of skipping
                // if we can not get a meaning value;
                if (j == i){
                    p_y_x[j] = -1;
                    continue;
                }
                
                // calculate p_y_x by filling temp_vec_var and calculate p_x by filling temp_vec_read
                int n_y_x = 0; int n_x = 0;
                for (int k = 0; k < pu_var[cand_loci[j]].size(); k++){
                    if (temp_vec_var[ pu_var[cand_loci[j]][k] ] == y_locus)
                        ++n_y_x;
                    if (temp_vec_read[ pu_var[cand_loci[j]][k] ] == y_read_locus)
                        ++n_x;
                }
                
                if (n_x < min_reads){
                    p_y_x[j] = -1;
                    continue;
                }
                p_y_x[j] = double(n_y_x) / n_x;
            }
            
            // sort p_y_x in descending order, and get index idx_p_y_x, i.e. p_y_x[idx_p_y_x[0]] is the maximum, 
            // p_y_x[idx_p_y_x[1]] is the second maximum and so on.
            vector<int> idx_p_y_x = sort_order(p_y_x, true); 
            
            // calculate conditional probability 
            int n_y_xp = 0; int n_xp = 0;
            int depth = 0;
            for (int j = 0; j < idx_p_y_x.size(); j++){
                if (cand_loci[ idx_p_y_x[j] ] == y_locus)
                    continue;
                if (depth > max_depth) break;
                
                int cur_locus = cand_loci[idx_p_y_x[j]];
                
                // calculate n_y_xp and n_xp
                n_y_xp = 0; n_xp = 0;
                if (j == 0) {
                    for (int k = 0; k < pu_var[cur_locus].size(); k++){
                        if (temp_vec_var[ pu_var[cur_locus][k] ] == y_locus){
                            temp_vec_var_lock[ pu_var[cur_locus][k] ] = cur_locus;
                            ++n_y_xp;
                        }
                        if (temp_vec_read[ pu_var[cur_locus][k] ] == y_read_locus){ 
                            temp_vec_read_lock[ pu_var[cur_locus][k] ] = cur_locus;
                            ++n_xp;
                        }
                    }
                }else{
                    for (int k = 0; k < pu_var[cur_locus].size(); k++){
                        if (temp_vec_var[ pu_var[cur_locus][k] ] == y_locus &&
                            temp_vec_var_lock[ pu_var[cur_locus][k] ] == cand_loci[idx_p_y_x[j-1]]){
                            temp_vec_var_lock[ pu_var[cur_locus][k] ] = cur_locus;
                            ++n_y_xp;
                        }
                        if (temp_vec_read[ pu_var[cur_locus][k] ] == y_read_locus &&
                            temp_vec_read_lock[ pu_var[cur_locus][k] ] == cand_loci[idx_p_y_x[j-1]]){ 
                            temp_vec_read_lock[ pu_var[cur_locus][k] ] = cur_locus;
                            ++n_xp;
                        }
                    }
                }
                
                if (n_xp < min_reads) break;
                
                double p_y_xp = (double)n_y_xp / n_xp;
                if (p_y_xp < cur_rl.p_y_xp)
                    break;
                
                // record result
                cur_rl.link_loci.push_back(cur_locus);
                cur_rl.n_y_xp = n_y_xp;
                cur_rl.n_xp = n_xp;
                cur_rl.p_y_xp = p_y_xp;
                
                ++depth; 
            }
            
            cur_rl.p_y_xp = double(n_y_xp) / n_xp;
            
            rl[y_locus].push_back(cur_rl);
            
        }
        
    }
    
};





#endif /* dforest_h */
