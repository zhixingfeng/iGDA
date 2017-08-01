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
    cout << "pileup reads" << endl;
    pu_var = pileup_var(encode_file, nreads_pu_var);
    pu_read = pileup_reads(align_file, nreads_pu_read);
    if (nreads_pu_var != nreads_pu_read)
        throw runtime_error("encode_file and align_file don't match.");
    
    cout << "subspace clustering" << endl;
    nreads = nreads_pu_var;
    // run the subpace clustering algorithm
    if (n_thread == 1){
        // single thread
        run_thread(cmpreads_file, out_file, max_cand_size, min_ratio, min_cvg, min_cvg);
    }else{
        // multiple threads
        // split the cmpreads_file
        cout << "split subspace" << endl;
        string tmp_prefix = tmp_dir + "/cmpreads_file_part"; 
        cmpreads_split(cmpreads_file, tmp_prefix, n_thread);
        
        // run with multiple threads
        cout << "run threads" << endl;
        vector<thread> threads;
        for (int i=0; i<n_thread; i++){
            string tmp_cmpreads_file = tmp_prefix + "_" + to_string(i);
            string tmp_out_file = tmp_dir + "/tmp_out_" + to_string(i) + ".sclust";
            threads.push_back(thread(&SClust::run_thread, this, tmp_cmpreads_file, tmp_out_file,
                                     max_cand_size, min_ratio, min_cvg, min_cvg));
        }
        for (int i=0; i<n_thread; i++)
            threads[i].join();

        // combine results
        string cmd = "cat ";
        for (int i=0; i<n_thread; i++){
            cmd += tmp_dir + "/tmp_out_" + to_string(i) + ".sclust ";
        }
        cmd += "> " + out_file;
        cout << cmd << endl;
        system(cmd.c_str());

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
            cout << "poccessed # of subspaces : " << k << endl;

        // load candidate subset
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_cmpreads_file);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_cmpreads_file);
        if (feof(p_cmpreads_file))
            break;
        
        // count frequency of variant combinations
        if (cand_loci.size() <= max_cand_size){
            // count frequency of pattern
            unordered_set<uint32_t> pattern; 
            int32_t nreads_cover_all = 0;
            this->count_freq(pattern, nreads_cover_all, cand_loci, temp_id_var, temp_id_read, temp_count_var);
            
            // skip if only few reads covering cand_loci
            if (nreads_cover_all >= min_cvg){
                // test each pattern for significance
                vector<uint32_t> rl_pattern;
                vector<int> rl_count;
                vector<double> rl_ratio;
                vector<double> rl_logLR;
                this->test_pattern(pattern, nreads_cover_all, temp_count_var, min_ratio, min_count, 
                                   rl_pattern, rl_logLR, rl_ratio, rl_count);
            
                // print frequency of pattern
                mtx_sclust.lock();
                print_pattern(p_outfile, cand_loci, rl_pattern, rl_logLR, rl_ratio, rl_count, nreads_cover_all);
                mtx_sclust.unlock();
            }
            // clear temp_count_var
            unordered_set<uint32_t>::iterator it;
            for (it=pattern.begin(); it!=pattern.end(); ++it)
                temp_count_var[*it] = 0;            
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
                        ++temp_count_var[temp_id_var[pu_read[var_read_locus][j]]];
                        pattern.insert(temp_id_var[pu_read[var_read_locus][j]]);
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

void SClust::test_pattern(unordered_set<uint32_t> &pattern, int32_t nreads_cover_all, vector<int32_t> &temp_count_var,
                          int min_ratio, int min_count, vector<uint32_t> &rl_pattern, 
                          vector<double> &rl_logLR, vector<double> &rl_ratio, vector<int> &rl_count)
{
    for (auto it = pattern.begin(); it != pattern.end(); ++it){
        // ignore single variants
        if (bitcount(*it)<=1)
            continue;
        
        // get significance of variant combinations
        if (temp_count_var[*it] >= min_count){
            double cur_min_ratio = double(nreads_cover_all) / temp_count_var[*it];
            double cur_min_logLR = cal_logLR(temp_count_var[*it], 0, 0, nreads_cover_all);
            //bool is_conditioned = false;
            for (auto it2 = pattern.begin(); it2 != pattern.end(); ++it2){
                if ((*it & *it2) == *it2 && *it > *it2){
                    //is_conditioned = true;
                    // get ratio between joint probability and product of marginal probability 
                    double cur_ratio = double(nreads_cover_all * temp_count_var[*it]) /
                                    ( double(temp_count_var[*it2] + temp_count_var[*it]) * double(temp_count_var[*it-*it2] + temp_count_var[*it]));
                    if (cur_ratio < cur_min_ratio)
                        cur_min_ratio = cur_ratio;
                    
                    // get likelihood ratio
                    double cur_logLR = cal_logLR(temp_count_var[*it], temp_count_var[*it2], temp_count_var[*it-*it2], nreads_cover_all);
                    if (cur_logLR < cur_min_logLR)
                        cur_min_logLR = cur_logLR;
                }
            }
            /*if (!is_conditioned){
                cur_min_ratio = double(nreads_cover_all) / temp_count_var[*it];
                //cur_min_logLR = 
            }*/
            
            if (cur_min_ratio >= min_ratio){
                rl_pattern.push_back(*it);
                rl_ratio.push_back(cur_min_ratio);
                rl_logLR.push_back(cur_min_logLR);
                rl_count.push_back(temp_count_var[*it]);
            }
        }
    }
}

void SClust::print_pattern(FILE *p_outfile, const vector<int> &cand_loci, vector<uint32_t> &rl_pattern,
                   vector<double> &rl_logLR, vector<double> &rl_ratio, vector<int> &rl_count, int32_t nreads_cover_all)
{
    
    for (int i=0; i<(int)rl_pattern.size(); ++i){
        // decode rl_pattern[i]
        bitset<32> pattern_bit(rl_pattern[i]);
        for (int j=0; j<(int)cand_loci.size(); ++j)
            fprintf(p_outfile, "%d,", cand_loci[j]);
        fprintf(p_outfile, "\t");
        
        // print pattern length
        int rl_pattern_len = bitcount(rl_pattern[i]);
        fprintf(p_outfile, "%d\t", rl_pattern_len);
        
        // print pattern_bit
        if (rl_pattern_len > cand_loci.size())
            throw runtime_error("incorrect # of 1s in rl_pattern");
        for (int j=0; j<(int)cand_loci.size(); ++j){
            if (pattern_bit[j]==1)
                fprintf(p_outfile, "%d,", cand_loci[j]);
        }
        fprintf(p_outfile, "\t%lf\t%lf\t%d\t%d\n", rl_ratio[i], rl_logLR[i], rl_count[i], nreads_cover_all);
    
        //fprintf(p_outfile, "\t%u\t%lf\t%lf\t%d\t%d\n", rl_pattern[i], rl_ratio[i], rl_logLR[i],
        //                                            rl_count[i], nreads_cover_all);
    }
}


// evaluate detected pattern
void SClust::eval_pattern(string pattern_file, string true_snp_file, string out_file)
{
    // load true snp file
    vector<unordered_set<int> > true_snp_data;
    ifstream fs_infile;
    open_infile(fs_infile, true_snp_file);
    int max_code = 0;
    while(1){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof()) break;
        vector<int> buf_vec = split_int(buf, ',');
        unordered_set<int> buf_vec_set;
        for (int i=0; i<(int)buf_vec.size(); ++i){
            if (buf_vec[i] > max_code)
                max_code = buf_vec[i];
        }
        buf_vec_set.insert(buf_vec.begin(), buf_vec.end());
        true_snp_data.push_back(buf_vec_set);
    }
    fs_infile.close();
    
    // initialize template vector
    vector<int> temp_vec(max_code + 1, 0);
    
    // scan pattern_file
    ofstream fs_outfile;
    open_outfile(fs_outfile, out_file);
    open_infile(fs_infile, pattern_file);
    while(1){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof()) break;
        vector<int> buf_vec = split_int(buf, ',');
        
        bool is_true = true;
        for (int i=0; i<(int)buf_vec.size(); ++i){
            bool is_in = false;
            for (int j=0; j<(int)true_snp_data.size(); ++j){
                if (true_snp_data[j].find(buf_vec[i]) != true_snp_data[j].end())
                    is_in = true;
            }
            if (!is_in)
                is_true = false;
        }
        fs_outfile << buf << '\t' << is_true << endl;
    }
    fs_infile.close();
    fs_outfile.close();
    
}


void SClust::summary(string sclust_file, string out_file, int min_overlap, double min_logLR)
{    
    ifstream fs_infile;    
    // scan the sclust file to get maximal code
    int max_code = 0;
    open_infile(fs_infile, sclust_file);
    while(1){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof()) break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=7)
            throw runtime_error("incorrect format in pattern_file");
        
        vector<int> pattern = split_int(buf_vec[2], ',');
        for (int i=0; i<(int)pattern.size(); ++i)
                if (pattern[i] > max_code)
                    max_code = pattern[i];
    }
    fs_infile.close();
    
    // the resulting subspaces 
    vector<vector<int> > subspaces; 
    
    // template to compare subspace and pattern
    vector<bool> temp_overlap(max_code, false);  
    vector<bool> temp_overlap_2(max_code, false);
    
    // scan the sclust file again 
    int64_t n_lines = 0;
    open_infile(fs_infile, sclust_file);
    while(1){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof()) break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=7)
            throw runtime_error("incorrect format in pattern_file");
        
        ++n_lines;
        if (n_lines%100000 == 0)
            cout << "n_line: " << n_lines << endl;

        if (stod(buf_vec[4]) < min_logLR)
            continue;
        vector<int> pattern = split_int(buf_vec[2], ',');
        
        if (subspaces.size()==0){
            // if subspaces is empty then add pattern as a new subspace
            subspaces.push_back(pattern);
        }else{
            // fill pattern to temp_overlap to compare it to subspaces
            for (int i=0; i<(int)pattern.size(); ++i)
                temp_overlap[pattern[i]] = true;
            
            // scan subspaces 
            bool is_exist = false;
            for (int i=0; i<(int)subspaces.size(); ++i){
                                    
                // count number of overlap
                int n_overlap = 0;
                for (int j=0; j<(int)subspaces[i].size(); ++j){
                    if (temp_overlap[subspaces[i][j]])
                        ++n_overlap;
                    temp_overlap_2[subspaces[i][j]] = true;
                }
                    
                // merge pattern and current subspace
                if (n_overlap >= min_overlap){
                    for (int j=0; j<(int)pattern.size(); ++j){
                        if (!temp_overlap_2[pattern[j]])
                            subspaces[i].push_back(pattern[j]);
                    }
                    is_exist = true;
                }
                
                // clear temp_overlap_2
                for (int j=0; j<(int)subspaces[i].size(); ++j)
                    temp_overlap_2[subspaces[i][j]] = false;
                
            }
            
            if (!is_exist)
                subspaces.push_back(pattern);
            
            // clear temp_overlap
            for (int i=0; i<(int)pattern.size(); ++i)
                temp_overlap[pattern[i]] = false;
        }
    }
    fs_infile.close();
    cout << "n_line: " << n_lines << endl;
    
    // output 
    ofstream fs_outfile;
    open_outfile(fs_outfile, out_file);
    for (int i=0; i<(int)subspaces.size(); ++i){
        for (int j=0; j<(int)subspaces[i].size(); ++j)
            fs_outfile << subspaces[i][j] << ",";
        fs_outfile << endl;
    }
    fs_outfile.close();
    
}



