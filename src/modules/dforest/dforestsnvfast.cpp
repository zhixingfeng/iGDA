//
//  dforestsnvfast.cpp
//  iGDA
//
//  Created by Zhixing Feng on 6/21/17.
//  Copyright (c) 2017 Zhixing Feng. All rights reserved.
//

#include "dforestsnvfast.h"

bool DForestSNVFast::run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread, double minfreq)
{
    cout << "number of threads: " << n_thread << endl;
    
    // load encode and alignment files
    cout << "pileup encode_file" << endl;
    call_pileup_var(encode_file);
    cout << "pileup align_file" << endl;
    call_pileup_reads(align_file);
    
    // initialize cache_n_y_x and cache_n_x
    for (int i=0; i<pu_var.size(); i++){
        cache_n_y_x.push_back(vector<int>(pu_var.size(), -1));
        cache_n_x.push_back(vector<int>(pu_var.size(), -1));
    }
    
    // single thread
    if (n_thread==1){
        run_thread(cmpreads_file, out_file, min_reads, max_depth, minfreq);
        return true;
    }

    return true;
}

void DForestSNVFast::build_tree(FILE * p_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq)
{
    // each of the locus in cand_loci is used as response y
    vector<double> p_y_x(cand_loci.size(), -1);
    for (int i = 0; i < cand_loci.size(); i++){
        DforestResult cur_rl;
        // NOTE: cand_loci[i] is response y. Let's fill in temp_vec_var and temp_vec_read
        // by reponse y. pu_var is 4 times larger than pu_read because of binary coding so
        // we have to devide 4 to access pu_read
        int y_locus = cand_loci[i];
        int y_read_locus = int (y_locus / 4);
        
        // fill in temp_vec_var and temp_vec_read by response y
        for (int j = 0; j < pu_var[y_locus].size(); j++)
            temp_vec_var[pu_var[y_locus][j]] = counter;
        
        for (int j = 0; j < pu_read[y_read_locus].size(); j++)
            temp_vec_read[pu_read[y_read_locus][j]] = counter;
        ++counter;
        if (counter >= numeric_limits<int64_t>::max()-1)
            throw runtime_error("counter exceeds maximal int64_t");
        
        // calculate joint frequency of neighbor variants
        for (int j = 0; j < cand_loci.size(); j++){
            // avoid self comparison. p_y_x will be reused, so give it -1 instead of skipping
            // if we can not get a meaning value;
            if (j == i){
                p_y_x[j] = -1;
                continue;
            }
            // check if n_y_x and n_x are cached 
            int n_y_x = 0; int n_x = 0;
            if (cache_n_y_x[cand_loci[i]][cand_loci[j]] == -1){
                // calculate p_y_x by filling temp_vec_var and calculate p_x by filling temp_vec_read
                for (int k = 0; k < pu_var[cand_loci[j]].size(); k++){
                    if (temp_vec_var[ pu_var[cand_loci[j]][k] ] == counter - 1)
                        ++n_y_x;
                    if (temp_vec_read[ pu_var[cand_loci[j]][k] ] == counter - 1)
                        ++n_x;
                }
                cache_n_y_x[cand_loci[i]][cand_loci[j]] = n_y_x;
                cache_n_x[cand_loci[i]][cand_loci[j]] = n_x;
            }else{
                n_y_x = cache_n_y_x[cand_loci[i]][cand_loci[j]];
                n_x = cache_n_x[cand_loci[i]][cand_loci[j]];
            }
            
            if (n_x < min_reads){
                p_y_x[j] = -1;
                continue;
            }
            if (n_x == 0)
                p_y_x[j] = 0;
            else
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
            if (depth >= max_depth) break;
            
            int cur_locus = cand_loci[idx_p_y_x[j]];
            
            // calculate n_y_xp and n_xp
            n_y_xp = 0; n_xp = 0;
            for (int k = 0; k < pu_var[cur_locus].size(); k++){
                if (temp_vec_var[ pu_var[cur_locus][k] ] == counter - 1){
                    temp_vec_var[ pu_var[cur_locus][k] ] = counter;
                    ++n_y_xp;
                }
                if (temp_vec_read[ pu_var[cur_locus][k] ] == counter - 1){
                    temp_vec_read[ pu_var[cur_locus][k] ] = counter;
                    ++n_xp;
                }
            }
            ++counter;
            if (counter >= numeric_limits<int64_t>::max()-1)
                throw runtime_error("counter exceeds maximal int64_t");
            
            if (n_xp < min_reads) break;
            
            double p_y_xp = (double)n_y_xp / n_xp;
            if (p_y_xp <= cur_rl.p_y_xp)
                break;
            
            // record result
            cur_rl.link_loci.push_back(cur_locus);
            cur_rl.n_y_xp = n_y_xp;
            cur_rl.n_xp = n_xp;
            cur_rl.p_y_xp = p_y_xp;
            
            ++depth; 
        }
        
        // write results (unordered) to outfile
        if (cur_rl.link_loci.size() > 0 && cur_rl.p_y_xp >= minfreq){
            fprintf(p_outfile, "%d\t%lf\t%lf\t%d\t%d\t%d\t", y_locus, cur_rl.bf, cur_rl.p_y_xp, cur_rl.n_y_xp, cur_rl.n_xp, (int)cur_rl.link_loci.size());
            for (int j = 0; j < cur_rl.link_loci.size(); j++)
                fprintf(p_outfile, "%d,", cur_rl.link_loci[j]);
            fprintf(p_outfile, "\n");
        }
    }
}

bool DForestSNVFast::run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq)
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