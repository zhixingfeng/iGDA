//
//  sclust.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/7/13.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "sclust.h"
mutex mtx_sclust;
void SClust::run(string encode_file, string align_file, string cmpreads_file, 
                 string out_file, string tmp_dir, int max_cand_size, int min_ratio, 
                 int min_count, int min_cvg, int n_thread)
{
    // initialize bit shift vector 
    for (int i=0; i<max_cand_size; ++i){
        bit_shift.push_back(uint32_t(pow(2,i)));
    }

    // pileup reads 
    int64_t nreads_pu_var, nreads_pu_read;
    pu_var = pileup_var(encode_file, nreads_pu_var);
    pu_read = pileup_reads(align_file, nreads_pu_read);
    if (nreads_pu_var != nreads_pu_read)
        throw runtime_error("encode_file and align_file don't match.");
    
    nreads = nreads_pu_var;
    // run the subpace clustering algorithm
    if (n_thread == 1){
        // single thread
        run_thread(cmpreads_file, out_file, max_cand_size, min_ratio, min_cvg, min_cvg);
    }else{
        // multiple threads
        
    }
}

bool SClust::run_thread(string cmpreads_file, string out_file, int max_cand_size, 
                        int min_ratio, int min_count, int min_cvg)
{
    // initialize templates
    vector<int32_t> temp_id_var(this->nreads, 0);
    vector<int32_t> temp_id_read(this->nreads, 0);
    vector<int32_t> temp_count_var(uint32_t(pow(2,max_cand_size)), 0);
    
    // open files
    FILE * p_cmpreads_file = fopen(cmpreads_file.c_str(), "rb");
    if (p_cmpreads_file == NULL)
        throw runtime_error("DForestSNVMax::run(): fail to open cmpreads_file");
    
    FILE *p_outfile = fopen(out_file.c_str(), "w");
    if (p_outfile == NULL)
        throw runtime_error("unable to open out_file");

    // scan cmpreads_file
    int64_t k = 1;
    while(1){
        if (k%10000==0)
            cout << "poccessed # of candidates : " << k << endl;

        // load candidate subset
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_cmpreads_file);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_cmpreads_file);
        if (feof(p_cmpreads_file))
            break;
        
        // count frequency of variant combinations
        if (cand_loci.size() <= max_cand_size){
            // count 
            unordered_set<uint32_t> pattern; 
            int32_t nreads_cover_all = 0;
            this->count_freq(pattern, nreads_cover_all, cand_loci, temp_id_var, temp_id_read, temp_count_var);
            
            // print frequency of pattern
            print_freq(p_outfile, cand_loci, pattern, nreads_cover_all, temp_count_var);
            
            // clear temp_count_var
            unordered_set<uint32_t>::iterator it;
            for (it=pattern.begin(); it!=pattern.end(); ++it){
                temp_count_var[*it] = 0;
            }
            
            
        }
        k++;
    }
    cout << "poccessed # of candidates : " << k << endl;
    
    fclose(p_cmpreads_file);
    fclose(p_outfile);

    
    return true;
}

void SClust::count_freq(unordered_set<uint32_t> &pattern, int32_t &nreads_cover_all,
                        const vector<int> &cand_loci, vector<int32_t> &temp_id_var, 
                        vector<int32_t> &temp_id_read, vector<int32_t> &temp_count_var)
{
    // get pattern covering all the candidate loci
    uint32_t pattern_cover_all = 0;
    for (int i=0; i<(int)cand_loci.size(); ++i)    
        pattern_cover_all += bit_shift[i];
    
    // scan candidate loci to find patterns of each read
    // reverse scan candidate loci because of binary encoding of pattern
    for (int i=0; i<(int)cand_loci.size(); ++i){
        // pu_var is 4 times larger than pu_read because of binary coding so we have to devide 4 to access pu_read
        int var_locus = cand_loci[i];
        int var_read_locus = int (var_locus / 4);
        
        // scan each locus to get variants and reads covering it
        for (int j=0; j<(int)pu_var[var_locus].size(); ++j)
            temp_id_var[pu_var[var_locus][j]] += bit_shift[i];
        
        for (int j=0; j<(int)pu_read[var_read_locus].size(); ++j){
            temp_id_read[pu_read[var_read_locus][j]] += bit_shift[i];
            // if hit the last the candidate locus, start to fill  temp_count_var and temp_count_var to count pattern        }
            if (i==cand_loci.size()-1){
                // only consider reads covering all the candidate loci
                if (temp_id_read[pu_read[var_read_locus][j]] == pattern_cover_all){
                    if (temp_id_var[pu_read[var_read_locus][j]] != 0){
                        ++temp_count_var[temp_id_var[ pu_read[var_read_locus][j]]];
                        pattern.insert(temp_id_var[ pu_read[var_read_locus][j]]);
                    }
                    ++nreads_cover_all;
                }                
            }
        }
    }
    
    // clear temp_id_var and temp_id_read;
    for (int i=0; i<(int)cand_loci.size(); ++i){
        
        int var_locus = cand_loci[i];
        int var_read_locus = int (var_locus / 4);
        
        for (int j=0; j<(int)pu_var[var_locus].size(); ++j)
            temp_id_var[pu_var[var_locus][j]] = 0;
        
        for (int j=0; j<(int)pu_read[var_read_locus].size(); ++j)
            temp_id_read[pu_read[var_read_locus][j]] = 0;
    }

}

void SClust::print_freq(FILE *p_outfile, const vector<int> &cand_loci, unordered_set<uint32_t> &pattern,
                int32_t nreads_cover_all, vector<int32_t> &temp_count_var)
{
    unordered_set<uint32_t>::iterator it;
    for (it=pattern.begin(); it!=pattern.end(); ++it){
        for (int i=0; i<(int)cand_loci.size(); ++i)
            fprintf(p_outfile, "%d,", cand_loci[i]);
        fprintf(p_outfile, "\t%u\t%d\t%d\n", *it, temp_count_var[*it], nreads_cover_all);
    }
}

