//
//  dforestsnvfast.cpp
//  iGDA
//
//  Created by Zhixing Feng on 6/21/17.
//  Copyright (c) 2017 Zhixing Feng. All rights reserved.
//

#include "dforestsnvfast.h"

bool DForestSNVFAST::run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread, double minfreq)
{
    cout << "number of threads: " << n_thread << endl;
    
    // load encode and alignment files
    cout << "pileup encode_file" << endl;
    call_pileup_var(encode_file);
    cout << "pileup align_file" << endl;
    call_pileup_reads(align_file);
    
    // initialize cache_n_y_x and cache_n_x
    for (int i=0; i<pu_var.size(); i++){
        cache_n_y_x.push_back(vector<int>(pu_var.size(),0));
        cache_n_x.push_back(vector<int>(pu_var.size(),0));
    }
    
    // single thread
    if (n_thread==1){
        run_thread(cmpreads_file, out_file, min_reads, max_depth, minfreq);
        return true;
    }

    return true;
}

void DForestSNVFAST::build_tree(FILE * p_cmpreads_file, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq)
{
    
}

bool DForestSNVFAST::run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq)
{
    // prepare buff of results and template
    vector<int64_t> temp_vec_var(this->n_reads, -1);
    vector<int64_t> temp_vec_read(this->n_reads, -1);
    
    
    // open cmpreads_file for each candidate subset, and output file
    int64_t k = 1;
    FILE * p_cmpreads_file = fopen(cmpreads_file.c_str(), "rb");
    if (p_cmpreads_file == NULL)
        throw runtime_error("DForestSNV::run(): fail to open cmpreads_file");
    
    FILE *p_outfile = fopen(out_file.c_str(), "w");
    if (p_outfile == NULL)
        throw runtime_error("unable to open out_file");
    
    // set counter and scan the candidates
    int64_t counter = 0;
    while(1){
        if (k%10000==0)
            cout << "poccessed # of candidates : " << k << endl;
        //printf("poccessed # of candidates : %d\n", k);
        // load candidate subset
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_cmpreads_file);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_cmpreads_file);
        if (feof(p_cmpreads_file))
            break;
        
        // build tree
        this->build_tree(p_outfile, cand_loci, counter, temp_vec_var, temp_vec_read, min_reads, max_depth, minfreq);
        
        k++;
    }
    cout << "poccessed # of candidates : " << k << endl;

    fclose(p_cmpreads_file);
    fclose(p_outfile);

    return true;
}
