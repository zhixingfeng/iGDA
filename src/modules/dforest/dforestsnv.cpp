//
//  dforestsnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/12/6.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "dforestsnv.h"
bool DForestSNV::run(string align_file, string encode_file, string cmpreads_file, string a_out_file, int min_reads, int max_depth)
{
    cout << "load encode_file" << endl;
    call_pileup_var(encode_file);
    cout << "load align_file" << endl;
    call_pileup_reads(align_file);
    
    out_file = a_out_file;
    
    // prepare buff of results and template 
    vector<int> temp_vec_var(this->n_reads, -1);
    vector<int> temp_vec_var_lock(this->n_reads, -1);
    vector<int> temp_vec_read(this->n_reads, -1);
    vector<int> temp_vec_read_lock(this->n_reads, -1);

    // scan cmpreads_file for each candidate subset
    int k = 1;
    FILE * p_cmpreads_file = fopen(cmpreads_file.c_str(), "rb");
    if (p_cmpreads_file == NULL)
        throw runtime_error("DForestSNV::run(): fail to open cmpreads_file");

    p_outfile = fopen(out_file.c_str(), "w");
    if (p_outfile == NULL)
        throw runtime_error("unable to open out_file");
    while(1){
        if (k%10000==0)
            printf("poccessed # of candidates : %d\n", k);
        // load candidate subset
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_cmpreads_file);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_cmpreads_file);
        if (feof(p_cmpreads_file))
            break;
        
        // build tree 
        build_tree(cand_loci, temp_vec_var, temp_vec_var_lock, temp_vec_read, temp_vec_read_lock, min_reads, max_depth);
        
        k++;
    }
    printf("poccessed # of candidates : %d\n", k);
    fclose(p_cmpreads_file);
    fclose(p_outfile);
    return true;
}

void DForestSNV::build_tree(const vector<int> &cand_loci, vector<int> &temp_vec_var, vector<int> &temp_vec_var_lock, vector<int> &temp_vec_read, vector<int> &temp_vec_read_lock, int min_reads, int max_depth)
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
        
        // write results (unordered) to outfile
        fprintf(p_outfile, "%d\t%lf\t%lf\t%d\t%d\t%d\t", y_locus, cur_rl.bf, cur_rl.p_y_xp, cur_rl.n_y_xp, cur_rl.n_xp, (int)cur_rl.link_loci.size());
        for (int j = 0; j < cur_rl.link_loci.size(); j++)
            fprintf(p_outfile, "%d,", cur_rl.link_loci[j]);
        fprintf(p_outfile, "\n");
            
    }
    
}
