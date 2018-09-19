//
//  assemble.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef assemble_h
#define assemble_h

#include <stxxl.h>
#include "../../../include/headers.h"
#include "../../misc/misc.h"
#include "../aligncoder/aligncodersnv.h"
#include "../src/modules/alignment/alignment.h"
#include "../src/modules/dforest/dforestsnvmax.h"

struct CmpreadsDiff
{
    CmpreadsDiff(){}
    CmpreadsDiff(const vector<int> &cand_loci, const vector<int> &cand_loci_diff):
            cand_loci(cand_loci), cand_loci_diff(cand_loci_diff),
            condprob(cand_loci.size(), -1), condprob_diff(cand_loci_diff.size(), -1),
            logLR(cand_loci.size(), -1), logLR_diff(cand_loci_diff.size(), -1){}
    int start;
    int end;
    const vector<int> cand_loci;
    const vector<int> cand_loci_diff;
    vector<double> condprob;
    vector<double> condprob_diff;
    vector<int> n_y_xp;
    vector<int> n_y_xp_diff;
    vector<double> logLR;
    vector<double> logLR_diff;
};
struct CmpreadsDiffRead
{
    CmpreadsDiffRead(){}
    CmpreadsDiffRead(int read_id): read_id(read_id){}
    int read_id;
    vector<CmpreadsDiff> cmpreads_diff;
    set<int> encode_corrected;
};

struct AdjEdge
{
    AdjEdge(){}
    AdjEdge(int left_node, int right_node, int n_match, int n_diff, double similarity): left_node(left_node), right_node(right_node),
            n_match(n_match), n_diff(n_diff), similarity(similarity){}
    int left_node;
    int right_node;
    int n_match;
    int n_diff;
    double similarity;
};


struct reads_compare
{
    bool operator()(const pair<int,double>& l, const pair<int,double> &r)
    {
        if (l.second != r.second){
            return l.second > r.second;
        }else{
            return l.first > r.first;
        }
    }
};


class Assembler
{
public:
    Assembler():cand_size(5),resampling_size(20),min_count(10),min_condprob(0.15),max_condprob(0.75){}
    Assembler(int cand_size, int resampling_size, int min_count, double min_condprob, double max_condprob):
            cand_size(cand_size),resampling_size(resampling_size),min_count(min_count),
            min_condprob(min_condprob),max_condprob(max_condprob){}
    virtual ~Assembler(){}
    
public:
    void get_variants(string dforest_file, string out_file, double min_condprob);
    void reduce_dim(string encode_file, string var_file, string out_file);
    void dist(string encode_file, string align_file, string out_file);
    void dist_rdim(string encode_file, string align_file, string var_file, string out_file);
    void jaccard_index(string encode_file, string align_file, string out_file, double min_jaccard_index);
    void jaccard_index_min(string encode_file, string align_file, string out_file, double cutoff);
    
    inline void call_pileup_var(string encode_file){
        pu_var.clear();
        pu_var = pileup_var(encode_file, n_reads);
    }
    
    inline void call_pileup_reads(string align_file, char format = 'm'){
        pu_read.clear();
        int64_t cur_n_reads;
        pu_read = pileup_reads(align_file, cur_n_reads, format);
        if (cur_n_reads != n_reads)
            throw runtime_error("number of reads in align_file and encode_file are different");
    }
    
    inline void set_par(int cand_size, int resampling_size, int min_count, double min_condprob, double max_condprob)
    {
        cand_size = cand_size;
        resampling_size = resampling_size;
        min_count = min_count;
        min_condprob = min_condprob;
        max_condprob = max_condprob;
    }

    // reconstruct reference genome from alignment
    void ref_reconstruct(const stxxl::vector<Align> &align_data, string &ref_name, string &ref_seq);
    
    // construct haplotype sequence
    void haplo_seq_construct(const vector<int> centroid, const string &ref_seq, string &haplo_seq);
    
    
    void run(string encode_file, string align_file, string out_file);
    
       
    /*----------- adaptive nearest neighbor clustering ------------*/
    // ann main function
    void ann_clust(string encode_file, string align_file, string var_file, int min_cvg = 20, double min_prop = 0.2, double max_prop = 0.7, int topn = 30, int max_nn = 200, double max_dist = 0.02);
    void ann_clust_recode(string recode_file, string recode_ref_file, string align_file, string var_file, int min_cvg = 20, double min_prop = 0.2, double max_prop = 0.7, int topn = 30, int max_nn = 200, double max_dist = 0.02);
    void find_nccontigs(vector<int64_t> &idx);
    void print_rl_ann_clust(string outfile, bool is_metric = false, vector<int64_t> idx = vector<int64_t>());
    void print_nc_reads_id(string outfile);
    
protected:
    // ann sub functions
    void find_ncreads(string encode_file, string align_file, string var_file, int topn = 30, double max_dist = 0.02);
    bool check_pileup(const vector<int> &pu_var_count, const vector<int> &pu_reads_count, int start, int end, const vector<int> &idx = vector<int>(), int min_cvg = 20, double min_prop = 0.2, double max_prop = 0.7);
    bool check_pileup_recode(const vector<int> &pu_var_count, const vector<int> &pu_var_ref_count, int start, int end, const vector<int> &idx = vector<int>(), int min_cvg = 20, double min_prop = 0.2, double max_prop = 0.7);
    void get_consensus_recode(ConsensusSeq &cons, const vector<int> &pu_var_count, const vector<int> &pu_var_ref_count, int start, int end, int min_cvg = 20);
    
    //void correct_contigs(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range, const vector<int> &var_cdf, const vector<bool> &temp_array, int min_cvg = 20, double min_prop = 0.2, double max_prop = 0.7);
    
    // output correct reads
    void print_correct_reads_raw(const CmpreadsDiffRead &cmpread, ofstream &fs_testfile);
    void print_correct_reads(const CmpreadsDiffRead &cmpread, ofstream &fs_outfile);
    
    
    
protected:
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_read;
    int64_t n_reads;
    
    int cand_size;
    int resampling_size;
    int min_count;
    double min_condprob;
    double max_condprob;
    
    
    vector<ConsensusSeq> rl_ann_clust;
    vector<int> nc_reads_id;
}
;


#endif /* assemble_h */


