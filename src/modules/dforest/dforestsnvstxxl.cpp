//
//  dforestsnvstxxl.cpp
//  iGDA
//
//  Created by Zhixing Feng on 2018/5/26.
//  Copyright © 2018年 Zhixing Feng. All rights reserved.
//

#include "dforestsnvstxxl.h"

bool DForestSNVSTXXL::run(const vector<vector<int> > &encode_data, const stxxl::vector<Align> &align_data,
         const stxxl::vector<vector<int> > &cmpreads_data, int min_reads, int max_depth,
         int n_thread, double minfreq)
{
    return true;
}

bool DForestSNVSTXXL::run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread, double minfreq, bool isinter)
{
    this->result.clear();
    this->result_all.clear();
    cout << "number of threads: " << n_thread << endl;
    
    // load encode and alignment files
    cout << "pileup encode_file" << endl;
    call_pileup_var(encode_file);
    
    cout << "pileup align_file" << endl;
    call_pileup_reads(align_file);
    
    cout << "filter encode_data" << endl;
    this->pu_var = filter_pileup_var(this->pu_var, this->pu_read, this->n_reads);
    
    // single thread
    if (n_thread==1){
        run_thread(cmpreads_file, out_file, min_reads, max_depth, minfreq, isinter);
    }else {
       
    }
    
    // save the result
    save_result(out_file, minfreq);
    //if (isinter)
        //save_result_all(out_file + ".all", minfreq);
    return true;
    
}

bool DForestSNVSTXXL::run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq, bool isinter)
{
    // prepare buff of results and template
    vector<int64_t> temp_vec_var(this->n_reads, -1);
    vector<int64_t> temp_vec_read(this->n_reads, -1);

    cout << "load cmpreads_data" << endl;
    stxxl::vector<vector<int> > cmpreads_data;
    loadcmpreads(cmpreads_data, cmpreads_file);
    
    cout << "build index for cmpreads_data" << endl;
    stxxl::vector<vector<int64_t> > cmpreads_index;
    build_index(cmpreads_index, cmpreads_data);
   
    cout << "build trees" << endl;
    ofstream fs_outfile;
    if (isinter)
        open_outfile(fs_outfile, out_file + ".all");
    int64_t counter = 0;
    result.resize(cmpreads_index.size());
    p_y_x_archive = vector<double>(cmpreads_index.size(),-1);
    int step_size = ceil(double(cmpreads_index.size())/100);
    for (auto i = 0; i < cmpreads_index.size(); ++i){
        if ((i+1)%step_size == 0)
            cout << "finished "  << floor(100*double(i)/cmpreads_index.size()) << "% = " << i << "/" << cmpreads_index.size() << endl;
        // build tree
        this->focal_locus = i;
        for (auto j = 0; j < cmpreads_index[i].size(); ++j){
            this->build_tree(fs_outfile, cmpreads_data[cmpreads_index[i][j]], counter, temp_vec_var, temp_vec_read, min_reads, max_depth, minfreq, isinter);
        }
        
        // clear p_y_x_archive
        for (auto it = idx_mod.begin(); it!=idx_mod.end(); ++it)
            p_y_x_archive[*it] = -1;
        idx_mod.clear();
    }
    cout << "finished 100% = " << cmpreads_index.size() << "/" << cmpreads_index.size() << endl;
    
    if (isinter)
        fs_outfile.close();
    return true;
}

void DForestSNVSTXXL::build_tree(ofstream &fs_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq, bool isinter)
{
    vector<double> p_y_x(cand_loci.size(), -1);
    DforestResult cur_rl;
    
    int y_locus = this->focal_locus;
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
        // if we can not get a meaningful value;
        if (cand_loci[j] == y_locus){
            p_y_x[j] = -1;
            continue;
        }
        
        if (p_y_x_archive[cand_loci[j]] != -1){
            p_y_x[j] = p_y_x_archive[cand_loci[j]];
        }else{
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
            
            p_y_x_archive[cand_loci[j]] = p_y_x[j];
            idx_mod.insert(cand_loci[j]);
        }
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
        
        // record current result
        cur_rl.focal_locus = y_locus;
        cur_rl.link_loci.push_back(cur_locus);
        cur_rl.n_y_xp = n_y_xp;
        cur_rl.n_xp = n_xp;
        cur_rl.p_y_xp = p_y_xp;
        
        ++depth;
    }
    
    // record result
    if (cur_rl.p_y_xp >= minfreq){
        // record maximal conditional probability of each locus
        if (cur_rl.p_y_xp > result[cur_rl.focal_locus].p_y_xp){
            result[cur_rl.focal_locus] = cur_rl;
        }else{
            // if p_y_xp equals, then prefer larger n_xp
            if (cur_rl.p_y_xp == result[cur_rl.focal_locus].p_y_xp){
                if (cur_rl.n_xp > result[cur_rl.focal_locus].n_xp){
                    result[cur_rl.focal_locus] = cur_rl;
                }else{
                    // if p_y_xp and n_xp equal, then prefer smaller link_loci size
                    if (cur_rl.n_xp == result[cur_rl.focal_locus].n_xp){
                        if (cur_rl.link_loci.size() < result[cur_rl.focal_locus].link_loci.size()){
                            result[cur_rl.focal_locus] = cur_rl;
                        }
                    }
                }
            }
        }
        
        // record result_all if isinter == true
        if (isinter){
            //result_all.push_back(cur_rl);
            if (cur_rl.link_loci.size() > 0){
                fs_outfile << cur_rl.focal_locus << '\t' << cur_rl.bf << '\t'
                << cur_rl.p_y_xp << '\t' << cur_rl.n_y_xp << '\t'
                << cur_rl.n_xp << '\t' << cur_rl.link_loci.size() << '\t'
                << cur_rl.link_loci << ',' << endl;
            }
        }
    }
}

void DForestSNVSTXXL::save_result(string out_file, double minfreq)
{
    ofstream fs_outfile;  open_outfile(fs_outfile, out_file);
    for (auto i = 0; i < result.size(); ++i){
        if (result[i].link_loci.size() > 0 && result[i].p_y_xp >= minfreq){
            fs_outfile << result[i].focal_locus << '\t' << result[i].bf << '\t'
            << result[i].p_y_xp << '\t' << result[i].n_y_xp << '\t'
            << result[i].n_xp << '\t' << result[i].link_loci.size() << '\t';
            for (auto j = 0; j < result[i].link_loci.size(); ++j)
                fs_outfile << result[i].link_loci[j] << ',';
            fs_outfile << endl;
        }
    }
    fs_outfile.close();
}

void DForestSNVSTXXL::save_result_all(string out_file, double minfreq)
{
    ofstream fs_outfile;  open_outfile(fs_outfile, out_file);
    for (auto i = 0; i < result_all.size(); ++i){
        if (result_all[i].link_loci.size() > 0 && result_all[i].p_y_xp >= minfreq){
            fs_outfile << result_all[i].focal_locus << '\t' << result_all[i].bf << '\t'
            << result_all[i].p_y_xp << '\t' << result_all[i].n_y_xp << '\t'
            << result_all[i].n_xp << '\t' << result_all[i].link_loci.size() << '\t';
            for (auto j = 0; j < result_all[i].link_loci.size(); ++j)
                fs_outfile << result_all[i].link_loci[j] << ',';
            fs_outfile << endl;
        }
    }
    fs_outfile.close();
}


void DForestSNVSTXXL::build_index(stxxl::vector<vector<int64_t> > &cmpreads_index, const stxxl::vector<vector<int> > &cmpreads_data)
{
    // get size of cmpreads_index
    int index_size = 0;
    for (auto i = 0; i < cmpreads_data.size(); ++i){
        for (auto j = 0; j < cmpreads_data[i].size(); ++j){
            if (cmpreads_data[i][j] + 1 >= index_size){
                index_size = cmpreads_data[i][j] + 1;
            }
        }
    }
    // build index
    cout << "index_size = " << index_size << endl;
    //for (int64_t i = 0; i < index_size; ++i)
    //    cmpreads_index.push_back(vector<int64_t>());
    //cmpreads_index.resize(index_size);
    cmpreads_index = stxxl::vector<vector<int64_t> >(index_size, 1);
    
    for (auto i = 0; i < cmpreads_data.size(); ++i){
        //if (i % 1000 == 0)
        //  cout << i << endl;
        for (auto j = 0; j < cmpreads_data[i].size(); ++j)
            cmpreads_index[cmpreads_data[i][j]].push_back(i);
    }
    
}




