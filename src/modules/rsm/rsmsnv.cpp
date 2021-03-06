//
//  rsmsnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 7/22/19.
//  Copyright © 2019 Zhixing Feng. All rights reserved.
//

#include "rsmsnv.h"


bool RSMsnv::run(string encode_file, string encode_ref_file, string align_file, string cmpreads_file, string out_file, int min_reads, int max_depth, int n_thread, double minfreq, double maxfreq, int min_homo_block_dist, bool isinter)
{
    cout << "number of threads: " << n_thread << endl;
    
    cout << "clear result" << endl;
    this->result.clear();
    
    cout << "pileup encode_file" << endl;
    int64_t n_reads_encode;
    this->pu_var.clear();
    this->pu_var = pileup_var(encode_file, n_reads_encode);
    
    cout << "pileup encode_ref_file" << endl;
    int64_t n_reads_encode_ref;
    this->pu_var_ref.clear();
    this->pu_var_ref = pileup_var(encode_ref_file, n_reads_encode_ref);
    
    if (n_reads_encode_ref != n_reads_encode)
        throw runtime_error("RSMsnv::run() incompactible encode_file and encode_ref_file");
    this->n_reads = n_reads_encode;
    
    cout << "construct pu_read" << endl;
    size_t pu_read_size = size_t((pu_var.size() - 1) / 4) + 1;
    this->pu_read.resize(pu_read_size);
    for (auto i = 0; i < pu_read.size(); ++i){
        for (int k = 0; k <= 3; ++k){
            if (4*i + k < pu_var.size())
                pu_read[i].insert(pu_read[i].end(), pu_var[4*i + k].begin(), pu_var[4*i + k].end());
            
            if (4*i + k < pu_var_ref.size())
                pu_read[i].insert(pu_read[i].end(), pu_var_ref[4*i + k].begin(), pu_var_ref[4*i + k].end());
        }
    }
    
    cout << "check homo_blocks" << endl;
    this->min_homo_block_dist = min_homo_block_dist;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    int64_t max_range = -1;
    for (auto i = 0; i < reads_range.size(); ++i)
        if (reads_range[i].second > max_range)
            max_range = reads_range[i].second;
    if (max_range + 1 > this->homo_blocks.size())
        throw runtime_error("Assembler::test_contigs, max_range + 1 > this->homo_blocks.size().");
    
    cout << "run random subspace maximization" << endl;
    if (n_thread==1){
        run_thread(cmpreads_file, out_file, min_reads, max_depth, minfreq, maxfreq, min_homo_block_dist, isinter);
    }else {
        
    }
    
    cout << "save results" << endl;
    save_result(out_file, minfreq);
    
    return true;
}

void RSMsnv::build_tree(ofstream &fs_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq, bool isinter)
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
        
        // exclude neighbors loci that are too close
        if (abs ( this->homo_blocks[int(cand_loci[j] / 4)] - this->homo_blocks[y_read_locus] ) < this->min_homo_block_dist){
            p_y_x[j] = -1;
            continue;
        }
        
        /*if (abs( int(cand_loci[j] / 4) - y_read_locus ) < 10){
         p_y_x[j] = -1;
         continue;
         }*/
        
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
    bool is_calculated = false;
    for (int j = 0; j < idx_p_y_x.size(); j++){
        if (cand_loci[ idx_p_y_x[j] ] == y_locus)
            continue;
        if (depth >= max_depth) break;
        if (p_y_x[idx_p_y_x[j]] < 0) break;
        
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
        
        is_calculated = true;
        ++depth;
    }
    
    // record result
    if (!is_calculated) return;
    
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
void RSMsnv::save_result(string out_file, double minfreq)
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


bool RSMsnv::run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq, double maxfreq, int min_homo_block_dist, bool isinter)
{
    // prepare buff of results and template
    vector<int64_t> temp_vec_var(this->n_reads, -1);
    vector<int64_t> temp_vec_read(this->n_reads, -1);
    
    cout << "load cmpreads_data" << endl;
    stxxl_vector_type_int cmpreads_data_raw;
    loadcmpreads(cmpreads_data_raw, cmpreads_file);
    const stxxl_vector_type_int &cmpreads_data = cmpreads_data_raw;
    
    cout << "build index for cmpreads_data" << endl;
    stxxl_vector_type cmpreads_index;
    build_index(cmpreads_index, cmpreads_data);
    
    size_t cmpreads_index_effective_size = 0;
    for (stxxl_vector_type::iterator it = cmpreads_index.begin(); it != cmpreads_index.end(); ++it)
        if (it->size() >0) ++cmpreads_index_effective_size;
    
    
    cout << "build trees" << endl;
    ofstream fs_outfile;
    if (isinter)
        open_outfile(fs_outfile, out_file + ".all");
    int64_t counter = 0;
    result.resize(cmpreads_index.size());
    p_y_x_archive = vector<double>(cmpreads_index.size(),-1);
    
    int step_size = ceil(double(cmpreads_index_effective_size)/100);
    
    int64_t i = 0;
    size_t n_scaned_loci = 0;
    for (stxxl_vector_type::iterator it = cmpreads_index.begin(); it != cmpreads_index.end(); ++it){
        if (it->size() > 0){
            ++n_scaned_loci;
            if (n_scaned_loci % step_size == 0)
                cout << "finished "  << floor(100*double(n_scaned_loci)/cmpreads_index_effective_size) << "% = " << n_scaned_loci << "/" << cmpreads_index_effective_size << endl;
        }
        // build tree
        this->focal_locus = (int)i;
        
        vector<int64_t> cur_cmpreads_index = *it;
        for (auto j = 0; j < cur_cmpreads_index.size(); ++j){
            this->build_tree(fs_outfile, cmpreads_data[cur_cmpreads_index[j]], counter, temp_vec_var, temp_vec_read, min_reads, max_depth, minfreq, isinter);
            // if p_y_xp is large enough then quit and test the next locus
            if (this->result[this->focal_locus].p_y_xp >= maxfreq)
                break;
        }
        
        // clear p_y_x_archive
        for (auto it = idx_mod.begin(); it!=idx_mod.end(); ++it)
            p_y_x_archive[*it] = -1;
        idx_mod.clear();
        
        ++i;
    }
    
    cout << "finished 100% = " << cmpreads_index.size() << "/" << cmpreads_index.size() << endl;
    
    if (isinter)
        fs_outfile.close();
    return true;
}

void RSMsnv::build_index(stxxl_vector_type &cmpreads_index, const stxxl_vector_type_int &cmpreads_data)
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
    
    cmpreads_index.resize(index_size);
    //cmpreads_index = cmpreads_index(index_size);
    
    for (auto i = 0; i < cmpreads_data.size(); ++i){
        //if (i % 1000 == 0)
        //  cout << i << endl;
        for (auto j = 0; j < cmpreads_data[i].size(); ++j)
            cmpreads_index[cmpreads_data[i][j]].push_back(i);
    }
}

