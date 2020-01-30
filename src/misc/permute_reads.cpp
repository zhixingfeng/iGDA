//
//  permute_reads.cpp
//  iGDA
//
//  Created by Zhixing Feng on 1/28/20.
//  Copyright © 2020 Zhixing Feng. All rights reserved.
//

#include "permute_reads.h"
#include "./io.h"

void permute_encodefile(string m5_file, string pileup_file, string outfile, int seed)
{
    // setup random number generator
    pcg32 rng(seed);
    
    // load range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, m5_file);
    
    // load pileup data
    unordered_map<int, vector<double> > pu_data;
    load_pileup(pu_data, pileup_file);
    
    // permute encode_data
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for (auto i = 0; i < reads_range.size(); ++i){
        for (auto locus = reads_range[i].first; locus <= reads_range[i].second; ++locus){
            // get substitution rate of the current locus
            vector<double> cur_freq = pu_data[locus];
            double sub_rate = 0;
            for (auto k = 0; k < cur_freq.size(); ++k)
                sub_rate += cur_freq[k];
            if (sub_rate < 0 || sub_rate > 1)
                throw runtime_error("permute_encodefile(): sub_rate < 0 || sub_rate > 1");
            
            // toss a coin to determine if there is a substitution
            vector<double> sub_prob = {1 - sub_rate, sub_rate};
            vector<int> sub_cand = {0, 1};
            int is_sub = gen_binom(sub_prob, sub_cand, rng);
            
            // randomly draw a substitution according to pileup data
            if (is_sub == 1){
                vector<int> cand = {0,1,2,3};
                int new_subtype = gen_binom(cur_freq, cand, rng);
                fs_outfile << 4*locus + new_subtype << "\t";
            }
        }
        fs_outfile << endl;
    }
    fs_outfile.close();
}

void get_condprob_threshold(string dforest_permuted_file, string pileup_file)
{
    unordered_map<int, double> dforest_data;
    load_dforestfile(dforest_data, dforest_permuted_file);
    
    unordered_map<int, vector<double> > pu_data;
    load_pileup(pu_data, pileup_file);
    
    int x = 1;
}
