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

bool DForestSNVSTXXL::run(string encode_file, string align_file, string cmpreads_file, string out_file, string tmp_dir, int min_reads, int max_depth, int n_thread, double minfreq)
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
       
    }
    
    // save the result
    save_result(out_file, minfreq);
    return true;
    
}

bool DForestSNVSTXXL::run_thread(string cmpreads_file, string out_file, int min_reads, int max_depth, double minfreq)
{
    cout << "load cmpreads_data" << endl;
    stxxl::vector<vector<int> > cmpreads_data;
    loadcmpreads(cmpreads_data, cmpreads_file);
    
    cout << "build index for cmpreads_data" << endl;
    vector<stxxl::vector<int64_t> > cmpreads_index;
    build_index(cmpreads_index, cmpreads_data);
    
    
    
    return true;
}

void DForestSNVSTXXL::build_tree(FILE * p_outfile, const vector<int> &cand_loci, int64_t &counter, vector<int64_t> &temp_vec_var, vector<int64_t> &temp_vec_read, int min_reads, int max_depth, double minfreq)
{
    
}


void DForestSNVSTXXL::save_result(string out_file, double minfreq)
{
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

}


void DForestSNVSTXXL::build_index(vector<stxxl::vector<int64_t> > &cmpreads_index, const stxxl::vector<vector<int> > &cmpreads_data)
{
    // get size of cmpreads_index
    size_t index_size = 0;
    for (auto i = 0; i < cmpreads_data.size(); ++i)
        for (auto j = 0; j < cmpreads_data[i].size(); ++j)
            if (cmpreads_data[i][j] + 1 >= index_size)
                index_size = cmpreads_data[i][j] + 1;

    // build index
    cmpreads_index.resize(index_size);
    for (auto i = 0; i < cmpreads_data.size(); ++i){
        for (auto j = 0; j < cmpreads_data[i].size(); ++j){
            cmpreads_index[cmpreads_data[i][j]].push_back(i);
        }
    }

}





