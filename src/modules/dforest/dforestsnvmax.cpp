//
//  dforestsnvmax.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/6/26.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "dforestsnvmax.h"
mutex mtx_snvmax;

/*---------input from memory/stxxl containers------------*/
bool DForestSNVMax::run(const vector<vector<int> > &encode_data, const stxxl::vector<Align> &align_data,
         const stxxl::vector<vector<int> > &cmpreads_data, int min_reads, int max_depth,
         int n_thread, double minfreq)
{
    this->result.clear();
    cout << "number of threads: " << n_thread << endl;
    
    // load encode and alignment files
    cout << "pileup encode_file" << endl;
    call_pileup_var(encode_data);
    cout << "pileup align_file" << endl;
    call_pileup_reads(align_data);

    // single thread
    if (n_thread==1){
        run_thread_stxxl(cmpreads_data, min_reads, max_depth, minfreq);
    }else {
        
    }
   
    return true;
}




/*-------------- input from files---------------*/
bool DForestSNVMax::run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread, double minfreq, bool isinter)
{
    this->result.clear();
    cout << "number of threads: " << n_thread << endl;
    
    // load encode and alignment files
    cout << "pileup encode_file" << endl;
    call_pileup_var(encode_file);
    cout << "pileup align_file" << endl;
    call_pileup_reads(align_file);
    
    // single thread
    if (n_thread==1){
        run_thread(cmpreads_file, out_file, min_reads, max_depth, minfreq);
    }else {
        /*// multiple threads
        // split the cmpreads_file
        string tmp_prefix = tmp_dir + "/cmpreads_file_part"; 
        cmpreads_split(cmpreads_file, tmp_prefix, n_thread);
    
        // run with multiple threads
        cout << "run threads" << endl;
        vector<thread> threads;
        for (int i=0; i<n_thread; i++){
            string tmp_cmpreads_file = tmp_prefix + "_" + to_string(i);
            string tmp_out_file = tmp_dir + "/tmp_out_" + to_string(i) + ".dforest";
            threads.push_back(thread(&DForestSNVMax::run_thread, this, tmp_cmpreads_file, tmp_out_file, min_reads, max_depth, minfreq));
        }
    
        for (int i=0; i<n_thread; i++)
            threads[i].join();*/
    }
    
    // write results (unordered) to outfile
    ofstream fs_outfile;  open_outfile(fs_outfile, out_file);
    for (it = result.begin(); it!=result.end(); ++it){
        if (it->second.link_loci.size() > 0 && it->second.p_y_xp >= minfreq){
            fs_outfile << it->second.focal_locus << '\t' << it->second.bf << '\t' 
                        << it->second.p_y_xp << '\t' << it->second.n_y_xp << '\t'
                        << it->second.n_xp << '\t' << it->second.link_loci.size() << '\t';
            for (int j = 0; j < it->second.link_loci.size(); j++)
                fs_outfile << it->second.link_loci[j] << ',';
            fs_outfile << endl;
        }
    }
    fs_outfile.close();
    return true;
    
}

void DForestSNVMax::build_tree(FILE * p_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq, bool isinter)
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
            
            // calculate p_y_x by filling temp_vec_var and calculate p_x by filling temp_vec_read
            int n_y_x = 0; int n_x = 0;
            for (int k = 0; k < pu_var[cand_loci[j]].size(); k++){
                if (temp_vec_var[ pu_var[cand_loci[j]][k] ] == counter - 1)
                    ++n_y_x;
                if (temp_vec_read[ pu_var[cand_loci[j]][k] ] == counter - 1)
                    ++n_x;
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
        mtx_snvmax.lock();
        vector<int> idx_p_y_x = sort_order(p_y_x, true); 
        mtx_snvmax.unlock();
        
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
            cur_rl.focal_locus = y_locus;
            cur_rl.link_loci.push_back(cur_locus);
            cur_rl.n_y_xp = n_y_xp;
            cur_rl.n_xp = n_xp;
            cur_rl.p_y_xp = p_y_xp;
            
            mtx_snvmax.lock();
            it = result.find(cur_rl.focal_locus);
            if (it != result.end()){
                // store cur_rl with larger p_y_xp 
                if (cur_rl.p_y_xp > it->second.p_y_xp){
                    it->second = cur_rl;
                }else{
                    if (cur_rl.p_y_xp == it->second.p_y_xp){
                        // if p_y_xp equals, then prefer larger n_xp
                        if (cur_rl.n_xp > it->second.n_xp){
                            it->second = cur_rl;
                        }else{
                            // if p_y_xp and n_xp equal, then prefer smaller link_loci size
                            if (cur_rl.n_xp == it->second.n_xp){
                                if (cur_rl.link_loci.size() < it->second.link_loci.size()){
                                    it->second = cur_rl;
                                }
                            }
                        }
                    }
                }
            }
            else{
                result[cur_rl.focal_locus] = cur_rl;
            }
            mtx_snvmax.unlock();
            ++depth; 
            
        }
        
    }
    
}

bool DForestSNVMax::run_thread_stxxl(const stxxl::vector<vector<int> > &cmpreads_data, int min_reads, int max_depth, double minfreq)
{
    // prepare buff of results and template
    vector<int64_t> temp_vec_var(this->n_reads, -1);
    vector<int64_t> temp_vec_read(this->n_reads, -1);
    
    // set counter and scan the candidates
    int64_t counter = 0;
    for (int64_t i=0; i<(int64_t)cmpreads_data.size(); ++i){
        if ((i+1) % 10000==0)
            cout << "poccessed # of candidates : " << i+1 << endl;
        //cout << cmpreads_data[i] << endl;
        this->build_tree(NULL, cmpreads_data[i], counter, temp_vec_var, temp_vec_read, min_reads, max_depth, minfreq);
    }
   
    cout << "poccessed # of candidates : " << cmpreads_data.size() << endl;
    
    
    return true;
}
bool DForestSNVMax::run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq, bool isinter)
{
    // prepare buff of results and template
    vector<int64_t> temp_vec_var(this->n_reads, -1);
    vector<int64_t> temp_vec_read(this->n_reads, -1);
    
    
    // open cmpreads_file for each candidate subset, and output file
    int64_t k = 1;
    FILE * p_cmpreads_file = fopen(cmpreads_file.c_str(), "rb");
    if (p_cmpreads_file == NULL)
        throw runtime_error("DForestSNVMax::run(): fail to open cmpreads_file");
    
    FILE *p_outfile = fopen(out_file.c_str(), "w");
    if (p_outfile == NULL)
        throw runtime_error("unable to open out_file");
    
    // set counter and scan the candidates    
    int64_t counter = 0;
    while(1){
        if (k%10000==0)
            cout << "poccessed # of candidates : " << k << endl;
        //printf("poccessed # of candidates : %d\n", (int)k);
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
