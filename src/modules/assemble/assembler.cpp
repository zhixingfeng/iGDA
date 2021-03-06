//
//  assemble.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "assembler.h"


void Assembler::get_variants(string dforest_file, string out_file, double min_condprob)
{
    AlignCoderSNV aligncodersnv;
    ifstream fs_dforeset_file;
    open_infile(fs_dforeset_file, dforest_file);
    ofstream fs_out_file;
    open_outfile(fs_out_file, out_file);
    while(1){
        string buf;
        getline(fs_dforeset_file, buf);
        if (fs_dforeset_file.eof())
            break;
        vector<string> buf_vec = split(buf,'\t');
        if (buf_vec.size()!=7)
            throw runtime_error("incorrect format in get_variants()");
        
        int code = stoi(buf_vec[0]);
        double condprob = stod(buf_vec[2]);
        if (condprob >= min_condprob){
            pair<int, char> rl_decode = aligncodersnv.decode(code);
            int locus = rl_decode.first;
            char base = rl_decode.second;
            fs_out_file << locus << '\t' << base << '\t' << buf << endl;
        }
        
    }
    fs_out_file.close();
    fs_dforeset_file.close();
    fs_out_file.close();
}

void Assembler::reduce_dim(string encode_file, string var_file, string out_file)
{
    // scan var_file to get maximal code
    ifstream fs_var_file;
    open_infile(fs_var_file, var_file);
    int max_code = 0;
    while (1){
        string buf;
        getline(fs_var_file, buf);
        if (fs_var_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in var_file");
        int cur_code = stoi(buf_vec[2]);
        if ( cur_code > max_code)
            max_code = cur_code;
    }
    fs_var_file.close();

    // fill template by scaning var_file for the second time
    vector<bool> temp_code(max_code + 1, false);
    open_infile(fs_var_file, var_file);
    while (1){
        string buf;
        getline(fs_var_file, buf);
        if (fs_var_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in var_file");
        int cur_code = stoi(buf_vec[2]);
        temp_code[cur_code] = true;
    }
    fs_var_file.close();
    
    // scan encode_file
    ifstream fs_encode_file;
    ofstream fs_out_file;
    open_infile(fs_encode_file, encode_file);
    open_outfile(fs_out_file, out_file);
    while (1){
        string buf;
        getline(fs_encode_file, buf);
        if (fs_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (int i=0; i<(int)buf_vec.size(); ++i){
            if (buf_vec[i] <= max_code){
                if (temp_code[buf_vec[i]]){
                    fs_out_file << buf_vec[i] << '\t';
                }
            }
        }
        fs_out_file << endl;
    }
    fs_var_file.close();
    fs_out_file.close();


}

void Assembler::dist(string encode_file, string align_file, string out_file)
{
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_file (m5 format)
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data and reads_range are different");
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0)
            cout << "processed " << i+1 << " reads" <<endl;
        
        // fill the template array by the variants in the ith read
        for (int k = 0; k < encode_data[i].size(); k++)
            temp_array[encode_data[i][k]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (i==j)
                continue;
            int overlap_start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            int overlap_end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = overlap_end - overlap_start + 1;
            //int n_overlap = (reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second) -
            //     (reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first) + 1;
            
            if (n_overlap <= 0)
                continue;
            
            // scan read i, only retain variants in [overlap_start, overlap_end]
            int n_miss = 0;
            for (int k = 0; k < encode_data[i].size(); k++){
                if (encode_data[i][k] >= overlap_start*4 &&
                    encode_data[i][k] <= overlap_end*4 + 3){
                    ++n_miss;
                }
            }
            
            // scan read j, only retain variants in [overlap_start, overlap_end]
            for (int k = 0; k < encode_data[j].size(); k++){
                if (encode_data[j][k] >= overlap_start*4 &&
                    encode_data[j][k] <= overlap_end*4 + 3){
                    if (temp_array[encode_data[j][k]] == i)
                        --n_miss;
                    else
                        ++n_miss;
                    
                }
            }
            
            // print results
            fs_out << i <<',' << j << ',' <<(double)n_miss / n_overlap << ',' << n_miss << ',' << n_overlap << endl;
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();
}


void Assembler::dist_rdim(string encode_file, string align_file, string var_file, string out_file)
{
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_file (m5 format)
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data and reads_range are different");
    
    // get maximal genome posotion
    int genome_size = 0;
    for (int i=0; i<(int)reads_range.size(); ++i)
        genome_size = genome_size < reads_range[i].second ? reads_range[i].second : genome_size;
    
    ++genome_size;
    
    // cumulated cdf of var
    ifstream fs_var_file;
    open_infile(fs_var_file, var_file);
    vector<int> cdf_var(genome_size, 0);
    int cur_count = 0;
    int pre_pos = 0;
    while (1){
        string buf;
        getline(fs_var_file, buf);
        if (fs_var_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in var_file");
        int cur_pos = stod(buf_vec[0]);
        
        for (int i=pre_pos; i<=cur_pos; ++i)
            cdf_var[i] = cur_count;
        
        ++cur_count;
        
        pre_pos = cur_pos + 1;
    }
    for (auto i = pre_pos; i < cdf_var.size(); ++i)
        cdf_var[i] = cur_count;
    fs_var_file.close();
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0)
            cout << "processed " << i+1 << " reads" <<endl;
        
        // fill the template array by the variants in the ith read
        for (int k = 0; k < encode_data[i].size(); k++)
            temp_array[encode_data[i][k]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (i==j)
                continue;
            int overlap_start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            int overlap_end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = overlap_end - overlap_start + 1;
           
            // calculate number of variants in [overlap_start, overlap_end]
            int n_var = cdf_var[overlap_end] - cdf_var[overlap_start];

            if (n_overlap <= 0 || cdf_var[overlap_end] - cdf_var[overlap_start] <=0)
                continue;
            
            // scan read i, only retain variants in [overlap_start, overlap_end]
            int n_miss = 0;
            for (int k = 0; k < encode_data[i].size(); k++){
                if (encode_data[i][k] >= overlap_start*4 &&
                    encode_data[i][k] <= overlap_end*4 + 3){
                    ++n_miss;
                }
            }
            
            // scan read j, only retain variants in [overlap_start, overlap_end]
            for (int k = 0; k < encode_data[j].size(); k++){
                if (encode_data[j][k] >= overlap_start*4 &&
                    encode_data[j][k] <= overlap_end*4 + 3){
                    if (temp_array[encode_data[j][k]] == i)
                        --n_miss;
                    else
                        ++n_miss;
                    
                }
            }
            
            // print results
            fs_out << i <<',' << j << ',' <<(double)n_miss / n_var << ',' << n_miss << ',' << n_overlap <<',' << n_var<< endl;
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();

}



void Assembler::jaccard_index(string encode_file, string align_file, string out_file, double min_jaccard_index)
{
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_file (m5 format)
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data and reads_range are different");
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size()-1; i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0)
            cout << "processed " << i+1 << " reads" <<endl;
        
        // fill the template array by the variants in the ith read
        for (int k = 0; k < encode_data[i].size(); k++)
            temp_array[encode_data[i][k]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=i+1; j<(int)encode_data.size(); j++){
            //if (i==j) continue;
            
            int start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            int end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = end - start + 1;
            
            int code_start = 4*start;
            int code_end = 4*end + 3;
            
            if (n_overlap <= 0)
                continue;
            
            int n_intersect = 0;
            int n_union = 0;
            for (int k = 0; k < encode_data[i].size(); k++){
                if (encode_data[i][k]>=code_start && encode_data[i][k]<=code_end)
                    ++n_union;
            }
            
            for (int k = 0; k < encode_data[j].size(); k++){
                if (encode_data[j][k] < code_start || encode_data[j][k] > code_end)
                    continue;
                if (temp_array[encode_data[j][k]] == i)
                    ++n_intersect;
                else
                    ++n_union;
            }
            
            double jaccard_index = n_union==0 ? 0 : (double)n_intersect / n_union;
            
            if (jaccard_index <= min_jaccard_index)
                continue;
            
            fs_out << i <<',' << j << ',' << jaccard_index << ',';
            fs_out << n_intersect << ',' << n_union << ',' << n_overlap << ',';
            fs_out << start << ',' << end << ',' << code_start << ',' << code_end << endl;
            
            /*fs_out << j <<',' << i << ',' << jaccard_index << ',';
            fs_out << n_intersect << ',' << n_union << ',' << n_overlap << ',';
            fs_out << start << ',' << end << ',' << code_start << ',' << code_end << endl;*/
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();
}

void Assembler::jaccard_index_min(string encode_file, string align_file, string out_file, double cutoff)
{
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_file (m5 format)
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file, 'm');
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data and reads_range are different");
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size()-1; i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0)
            cout << "processed " << i+1 << " reads" <<endl;
        
        // fill the template array by the variants in the ith read
        for (int k = 0; k < encode_data[i].size(); k++)
            temp_array[encode_data[i][k]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=i+1; j<(int)encode_data.size(); j++){
            if (reads_range[i].first > reads_range[j].second || reads_range[j].first > reads_range[i].second)
                continue;
            
            int n_intersect = 0;
            for (int k = 0; k < encode_data[j].size(); k++)
                if (temp_array[encode_data[j][k]] == i)
                    ++n_intersect;
            double condProb_i = encode_data[i].size() == 0 ? 0 : (double)n_intersect / encode_data[i].size();
            double condProb_j = encode_data[j].size() == 0 ? 0 : (double)n_intersect / encode_data[j].size();
            
            double condProb = condProb_i <= condProb_j ? condProb_i : condProb_j;
            
            if (condProb < cutoff)
                continue;
            
            fs_out << i <<',' << j << ',' << condProb << ',';
            fs_out << condProb_i << ',' << condProb_j << ',';
            fs_out << n_intersect << ',' << encode_data[i].size() << ',' << encode_data[j].size() << endl;
            
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();
}


void Assembler::ref_reconstruct(const vector<Align> &align_data, string &ref_name, string &ref_seq)
{
    // get reference genome name (only 1 chromosome is allowed)
    if (align_data.size() == 0)
        throw runtime_error("align_data is empty.");
    ref_name = align_data[0].tName;
    
    // scan align_data to get genome size and check if chromosome is unique
    int g_size = -1;
    for (int i=0; i<(int)align_data.size(); ++i){
        if (align_data[0].tName != ref_name)
            runtime_error("ref_name is not unique");
        if (align_data[i].tEnd > g_size)
            g_size = align_data[i].tEnd;
    }
    g_size = g_size + 1;
    
    // reconstruct ref_seq
    ref_seq = string ('N', g_size);
    for (int i=0; i<(int)align_data.size(); ++i){
        string cur_tSeq = align_data[i].tSeq;
        if (cur_tSeq.size() != align_data[i].tEnd - align_data[i].tStart + 1)
            throw runtime_error("cur_tSeq.size() != align_data[i].tEnd - align_data[i].tStart + 1");
        ref_seq.replace(align_data[i].tStart, cur_tSeq.size(), cur_tSeq);
    }
}

void Assembler::haplo_seq_construct(const vector<int> centroid, const string &ref_seq, string &haplo_seq)
{
    AlignCoderSNV aligncoder;
    haplo_seq = ref_seq;
    for (int i=0; i<(int)centroid.size(); ++i){
        pair<int, char> cur_decode = aligncoder.decode(centroid[i]);
        if (cur_decode.first >= haplo_seq.size())
            throw runtime_error("cur_decode.first >= haplo_seq.size()");
        haplo_seq[cur_decode.first] = cur_decode.second;
    }
}

void Assembler::run(string encode_file, string align_file, string out_file)
{
    

    
}


void Assembler::ann_clust(string encode_file, string align_file, string var_file, int min_cvg, double min_prop, double max_prop, int topn, int max_nn, double max_dist)
{
    /*------------ find nc-reads -----------*/
    cout << "find non-contained reads" << endl;
    this->find_ncreads(encode_file, align_file, var_file, topn, max_dist);
    cout << "number of nc-reads: " << nc_reads_id.size() << endl;
    
    /*------------ use nc-reads seed to cluster ----------*/
    cout << "use non-contained reads as seed to cluster" << endl;
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (reads_range.size() != encode_data.size())
        throw runtime_error("reads_range.size() != encode_data.size()");
    
    // to be removed
    //nc_reads_id.resize(encode_data.size());
    //iota(nc_reads_id.begin(), nc_reads_id.end(), 0);
    
    // get genome size
    size_t genome_size = get_genome_size(reads_range);
    
    // get cumulative sum of variants
    vector<int> var_cdf; get_var_cdf(var_cdf, var_file, genome_size);
    
    // create a template to compare reads
    vector<bool> temp_array(genome_size*4+3, false);
    
    vector<int> cur_pu_var_count(4*(genome_size-1)+3+1, 0);
    vector<int> cur_pu_reads_count(genome_size, 0);
    
    // get topn nearest neighbors for each reads and find non-contained reads (no topn reads can cover all its range)
    // to be removed
    //ofstream fs_outfile;
    //open_outfile(fs_outfile, encode_file + ".igda_dist");

    int64_t n_nc_reads = 0;
    for (auto i : nc_reads_id){
        ++n_nc_reads;
        if (n_nc_reads % 1000 == 0)
            cout << "processed " << n_nc_reads << " / " << nc_reads_id.size() << endl;
        // calculate hamming distance between reads i and other reads
        priority_queue<pair<int,double>, vector<pair<int,double> >, reads_compare_dist > topn_id;
        for (auto j = 0; j < encode_data.size(); ++j){
            if (i == j) continue;
            if (reads_range[i].first >= reads_range[j].second || reads_range[i].second <= reads_range[j].first)
                continue;
            if (reads_range[i].first < reads_range[j].first)
                continue;
            
            double cur_dist = dist_hamming(encode_data[i], encode_data[j], reads_range[i], reads_range[j], var_cdf, temp_array);
            
            if (cur_dist < 0) continue;
            
            topn_id.push(pair<int,double>(j,cur_dist));
            
            // to be removed
            //fs_outfile << i << '\t' << j << '\t' << cur_dist << endl;
        }

        // get max_nn neighbors
        vector<int> cur_neighbors_topn;
        vector<int> cur_neighbors;
        vector<double> cur_neighbors_dist;
        for (auto j = 0; j < max_nn; ++j){
            if (topn_id.empty()) break;
            int cur_id = topn_id.top().first;
            double cur_dist = topn_id.top().second;
            cur_neighbors.push_back(cur_id);
            cur_neighbors_dist.push_back(cur_dist);
            if (j < topn)
                cur_neighbors_topn.push_back(cur_id);
            topn_id.pop();
        }
        
        if (cur_neighbors.size() == 0 || cur_neighbors_topn.size() == 0)
            continue;
        
        // pileup all neighbors and pop from most distant neighbor until all loci are homogeneous
        // pileup topn neighbors
        
        unordered_set<int64_t> mod_idx_var;
        ReadRange mod_range(reads_range[cur_neighbors[0]]);
        for (auto j : cur_neighbors){
            pileup_var_online_count(cur_pu_var_count, encode_data[j], mod_idx_var);
            pileup_reads_m5_online_count(cur_pu_reads_count, reads_range[j], mod_range);
        }
        
        // check pileup of all the neighbors
        bool is_homo = this->check_pileup(cur_pu_var_count, cur_pu_reads_count, reads_range[i].first, reads_range[i].second, vector<int>(), min_cvg, min_prop, max_prop);
        
        // if not all loci are homogeneous, pop neighbors from the most distant one until all loci are homogeneous or number of neighbors <= topn
        if (!is_homo){
            for (auto j = 0; j < cur_neighbors.size(); ++j){
                int t = (int)cur_neighbors.size() - 1 - j;
                if (t < topn) break;
                pileup_var_online_count_pop(cur_pu_var_count, encode_data[cur_neighbors[t]]);
                pileup_reads_m5_online_count_pop(cur_pu_reads_count, reads_range[cur_neighbors[t]]);
                is_homo = this->check_pileup(cur_pu_var_count, cur_pu_reads_count, reads_range[i].first, reads_range[i].second, vector<int>(), min_cvg, min_prop, max_prop);
                if (is_homo) break;
            }
        }

        // if we can get all loci homogeneous get consensus sequence
        if (is_homo){
            ConsensusSeq cur_cons;
            cur_cons.seed = encode_data[i];
            cur_cons.neighbors_id = cur_neighbors;
            get_consensus(cur_cons, cur_pu_var_count, cur_pu_reads_count, reads_range[i].first, reads_range[i].second, min_cvg);
            rl_ann_clust.push_back(cur_cons);
        }
        
        // clean cur_pu_var_count and cur_pu_reads_count
        for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
            cur_pu_var_count[*it] = 0;
        for (auto j = mod_range.first; j <= mod_range.second; ++j)
            cur_pu_reads_count[j] = 0;
        
        
    }
    // to be removed
    //fs_outfile.close();
    //cout << "processed " << nc_reads_id.size() << " / " << nc_reads_id.size() << endl;
    
    // correct contigs
    //this->correct_contigs(encode_data, reads_range, var_cdf, temp_array, min_cvg, min_prop, max_prop);
    
}

void Assembler::ann_clust_recode(string recode_file, string recode_ref_file, string encode_file, string align_file, string var_file, int min_cvg, double min_prop, double max_prop, int topn, int max_nn, double min_jaccard, bool is_correct, bool is_hang, int max_iter, bool is_recode, vector<int64_t> reads_id)
{
    /*------------ find nc-reads (deperated) -----------*/
    //cout << "find non-contained reads" << endl;
    //this->find_ncreads(recode_file, align_file, var_file, topn, 0.02);
    //cout << "number of nc-reads: " << nc_reads_id.size() << endl;
    
    /*------------ use nc-reads seed to cluster ----------*/
    //cout << "use non-contained reads as seed to cluster" << endl;
    
    // load encode data
    cout << "load encode data" << endl;
    vector<vector<int> > encode_data;
    vector<vector<int> > encode_data_cd;
    loadencodedata(encode_data, encode_file);

    // load recode data
    cout << "load recode data" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    // load recode_ref data
    cout << "load recode_ref data" << endl;
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_ref_file);
    
    // load reads range
    cout << "load m5 data" << endl;
    vector<ReadRange> reads_range;
    vector<ReadRange> reads_range_cd;
    loadreadsrange(reads_range, align_file);
    
    if (reads_range.size() != encode_data.size())
        throw runtime_error("reads_range.size() != encode_data.size()");
    
    if (reads_range.size() != recode_data.size())
        throw runtime_error("reads_range.size() != recode_data.size()");
    
    if (reads_range.size() != recode_ref_data.size())
        throw runtime_error("reads_range.size() != recode_ref_data.size()");
    
    if (reads_id.size() == 0){
        this->nc_reads_id.resize(reads_range.size());
        iota(this->nc_reads_id.begin(), this->nc_reads_id.end(), 0);
    }else{
        this->nc_reads_id = reads_id;
    }
    
    // to be removed
    //this->nc_reads_id = {14831};
    
    cout << "number of nc-reads: " << nc_reads_id.size() << endl;
    
    
    // get genome size
    size_t genome_size = get_genome_size(reads_range);
    
    // get cumulative sum of variants
    //vector<int> var_cdf; get_var_cdf(var_cdf, var_file, genome_size);
    
    // create a template to compare reads
    vector<bool> temp_array(genome_size*4+3, false);
    
    vector<int> cur_pu_var_count(4*(genome_size-1)+3+1, 0);
    vector<int> cur_pu_var_ref_count(4*(genome_size-1)+3+1, 0);
    vector<int> cur_pu_reads_count(genome_size, 0);
    
    // get topn nearest neighbors for each reads and find non-contained reads (no topn reads can cover all its range)
    if (is_correct){
        cout << "correct reads" << endl;
        correct_reads(encode_data, reads_range, encode_data_cd, reads_range_cd);
    }
    ofstream fs_iterfile;
    open_outfile(fs_iterfile, encode_file + ".iter");
    int64_t n_nc_reads = 0;
    for (auto i : nc_reads_id){
        ++n_nc_reads;
        if (n_nc_reads % 1000 == 0)
            cout << "processed " << n_nc_reads << " / " << nc_reads_id.size() << endl;
        
        ConsensusSeq cur_cons;
        if (is_correct){
            cur_cons.seed = encode_data_cd[i];
            cur_cons.cons_seq = encode_data_cd[i];
            cur_cons.start = reads_range_cd[i].first;
            cur_cons.end = reads_range_cd[i].second;
        }else{
            if (is_recode){
                cur_cons.seed = recode_data[i];
                cur_cons.cons_seq = recode_data[i];
            }else{
                cur_cons.seed = encode_data[i];
                cur_cons.cons_seq = encode_data[i];
            }
            cur_cons.start = reads_range[i].first;
            cur_cons.end = reads_range[i].second;
        }
        bool is_valid = false;
        vector<int> prev_cons_seq = cur_cons.cons_seq;
        int b = 0;
        for (b = 0; b < max_iter; ++b){
            // calculate hamming distance between reads i and other reads and get topn nearest neighbors
            //priority_queue<pair<int,double>, vector<pair<int,double> >, reads_compare_dist > topn_id;
            priority_queue<pair<int,double>, vector<pair<int,double> >, reads_compare_sim > topn_id;
            for (auto j = 0; j < recode_data.size(); ++j){
                if (i == j) continue;
                
                if (cur_cons.start >= reads_range[j].second || cur_cons.end <= reads_range[j].first)
                    continue;
                
                if (cur_cons.start < reads_range[j].first && is_hang)
                    continue;
                
                //double cur_dist = dist_hamming(recode_data[i], recode_data[j], reads_range[i], reads_range[j], var_cdf, temp_array);
                //double cur_dist = sim_jaccard(cur_cons.cons_seq, recode_data[j], reads_range[i], reads_range[j], temp_array, true);
                int min_overlap = int(0.75*(cur_cons.end - cur_cons.start));
                double cur_dist;
                if (is_correct){
                    cur_dist = sim_jaccard(encode_data_cd[i], encode_data[j], reads_range_cd[i], reads_range[j], temp_array, true, min_overlap);
                }else{
                    ReadRange cur_range(cur_cons.start, cur_cons.end);
                    if (is_recode){
                        cur_dist = sim_jaccard(cur_cons.cons_seq, recode_data[j], cur_range, reads_range[j], temp_array, true, min_overlap, true);
                    }else{
                        cur_dist = sim_jaccard(cur_cons.cons_seq, encode_data[j], cur_range, reads_range[j], temp_array, true, min_overlap, true);
                    }
                    //cur_dist = sim_jaccard(cur_cons.cons_seq, recode_data[j], cur_range, reads_range[j], temp_array, true, 0, true);
                    //cur_dist = sim_jaccard(encode_data[i], encode_data[j], reads_range[i], reads_range[j], temp_array, true, 0, true);
                    //cur_dist = sim_jaccard(encode_data[i], encode_data[j], reads_range[i], reads_range[j], temp_array, true, min_overlap);
                }
                if (cur_dist <= min_jaccard) continue;
                
                topn_id.push(pair<int,double>(j,cur_dist));
            }
            
            // get max_nn neighbors
            //vector<int> cur_neighbors_topn;
            vector<int> cur_neighbors;
            vector<double> cur_neighbors_dist;
            for (auto j = 0; j < max_nn; ++j){
                if (topn_id.empty()) break;
                int cur_id = topn_id.top().first;
                double cur_dist = topn_id.top().second;
                cur_neighbors.push_back(cur_id);
                cur_neighbors_dist.push_back(cur_dist);
                //if (j < topn)
                //    cur_neighbors_topn.push_back(cur_id);
                topn_id.pop();
            }
            
            if (cur_neighbors.size() < topn || cur_neighbors.size() < min_cvg)
                break;
            
            // pileup all neighbors and pop from most distant neighbor until all loci are homogeneous
            // pileup topn neighbors
            
            unordered_set<int64_t> mod_idx_var;
            unordered_set<int64_t> mod_idx_var_ref;
            for (auto j : cur_neighbors){
                pileup_var_online_count(cur_pu_var_count, recode_data[j], mod_idx_var);
                pileup_var_online_count(cur_pu_var_ref_count, recode_ref_data[j], mod_idx_var_ref);
            }
            
            // check pileup of all the neighbors
            bool is_homo = this->check_pileup_recode(cur_pu_var_count, cur_pu_var_ref_count, cur_cons.start, cur_cons.end, vector<int>(), min_cvg, min_prop, max_prop);
            
            // if not all loci are homogeneous, pop neighbors from the most distant one until all loci are homogeneous or number of neighbors <= topn
            int t = (int)cur_neighbors.size();
            if (!is_homo){
                for (auto j = 0; j < cur_neighbors.size(); ++j){
                    if (t <= topn) break;
                    t = (int)cur_neighbors.size() - 1 - j;
                    pileup_var_online_count_pop(cur_pu_var_count, recode_data[cur_neighbors[t]]);
                    pileup_var_online_count_pop(cur_pu_var_ref_count, recode_ref_data[cur_neighbors[t]]);
                    is_homo = this->check_pileup_recode(cur_pu_var_count, cur_pu_var_ref_count, cur_cons.start, cur_cons.end, vector<int>(), min_cvg, min_prop, max_prop);
                    if (is_homo) break;
                }
            }
            cur_cons.neighbors_id = vector<int> (cur_neighbors.begin(),  cur_neighbors.begin() + t);
            
            // to be removed
            //this->get_consensus_recode(cur_cons, cur_pu_var_count, cur_pu_var_ref_count, cur_cons.start, cur_cons.end, min_cvg);
            //cout << "cur_cons.cons_seq: " << cur_cons.cons_seq << endl;
            //cout << "cur_cons.neighbors_id: " << cur_cons.neighbors_id << endl;
            
            if (is_homo){
                this->get_consensus_recode(cur_cons, cur_pu_var_count, cur_pu_var_ref_count, cur_cons.start, cur_cons.end, min_cvg);
                is_valid = true;
                
                if (std::equal(cur_cons.cons_seq.begin(), cur_cons.cons_seq.end(), prev_cons_seq.begin())){
                    // clean cur_pu_var_count and cur_pu_reads_count
                    for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
                        cur_pu_var_count[*it] = 0;
                    for (auto it = mod_idx_var_ref.begin(); it != mod_idx_var_ref.end(); ++it)
                        cur_pu_var_ref_count[*it] = 0;
                    
                    break;
                }else{
                    prev_cons_seq = cur_cons.cons_seq;
                }
            }else{
                is_valid = false;
                
                // clean cur_pu_var_count and cur_pu_reads_count
                for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
                    cur_pu_var_count[*it] = 0;
                for (auto it = mod_idx_var_ref.begin(); it != mod_idx_var_ref.end(); ++it)
                    cur_pu_var_ref_count[*it] = 0;
                
                break;
            }
            
            // clean cur_pu_var_count and cur_pu_reads_count
            for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
                cur_pu_var_count[*it] = 0;
            for (auto it = mod_idx_var_ref.begin(); it != mod_idx_var_ref.end(); ++it)
                cur_pu_var_ref_count[*it] = 0;
            
            //if (is_homo)
            //    break;
        }
        if (is_valid){
            rl_ann_clust.push_back(cur_cons);
            fs_iterfile << cur_cons.cons_seq << "\t" << b + 1 << endl;
        }
    }
    fs_iterfile.close();
    
}

void Assembler::ann_clust_recode_multithread(string recode_file, string recode_ref_file, string encode_file, string align_file,
                                  int min_cvg, double min_prop, double max_prop, int topn,
                                  int max_nn, double min_jaccard, bool is_correct, bool is_hang,
                                  int max_iter, bool is_recode, int nthread)
{
    cout << "load encode data" << endl;
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load recode data
    cout << "load recode data" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    // load recode_ref data
    cout << "load recode_ref data" << endl;
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_ref_file);
    
    // load reads range
    cout << "load m5 data" << endl;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (reads_range.size() != encode_data.size())
    throw runtime_error("reads_range.size() != encode_data.size()");
    
    if (reads_range.size() != recode_data.size())
    throw runtime_error("reads_range.size() != recode_data.size()");
    
    if (reads_range.size() != recode_ref_data.size())
    throw runtime_error("reads_range.size() != recode_ref_data.size()");
    
    // split reads and run multithread
    size_t blocksize_leftover = encode_data.size() % nthread;
    size_t blocksize = encode_data.size() / nthread;
    
    if (blocksize == 0)
    throw runtime_error("cmpreads_topn_multithread(): blocksize == 0");
    
    size_t cur_index = 0;
    vector<pair<size_t, size_t> > block_range;
    while(cur_index < encode_data.size()){
        if (cur_index + blocksize <= encode_data.size()){
            block_range.push_back(pair<size_t, size_t>(cur_index, cur_index + blocksize - 1));
        }else{
            if (blocksize_leftover == 0)
                throw runtime_error("cmpreads_topn_multithread(): blocksize_leftover == 0");
            block_range.back().second += blocksize_leftover;
        }
        cur_index += blocksize;
    }
    
    // run multithread
    vector<thread> threads;
    for (auto i = 0; i < block_range.size(); ++i){
        if (block_range[i].first < 0 || block_range[i].second >= encode_data.size())
            throw runtime_error("cmpreads_topn_multithread(): block_range[i].first < 0 || block_range[i].second >= encode_data.size()");
        vector<int64_t> cur_idx;
        for (auto j = block_range[i].first; j <= block_range[i].second; ++j)
            cur_idx.push_back(j);
        
        threads.push_back(thread(&Assembler::ann_clust_recode_core, this, std::ref(recode_data), std::ref(recode_ref_data), std::ref(encode_data), std::ref(reads_range), min_cvg,
                                 min_prop, max_prop, topn, max_nn, min_jaccard, is_correct, is_hang, max_iter, is_recode, cur_idx));
    }
    
    for (auto i = 0; i < threads.size(); ++i)
        threads[i].join();
    
    //int x = 1;
}

void Assembler::ann_clust_recode_core(const vector<vector<int> > &recode_data, const vector<vector<int> > &recode_ref_data,
                      const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                      int min_cvg, double min_prop, double max_prop, int topn,
                      int max_nn, double min_jaccard, bool is_correct, bool is_hang,
                      int max_iter, bool is_recode, vector<int64_t> reads_id)
{
    if (reads_range.size() != encode_data.size())
    throw runtime_error("reads_range.size() != encode_data.size()");
    
    if (reads_range.size() != recode_data.size())
    throw runtime_error("reads_range.size() != recode_data.size()");
    
    if (reads_range.size() != recode_ref_data.size())
    throw runtime_error("reads_range.size() != recode_ref_data.size()");
    
    //this->nc_reads_id = reads_id;
    cout << "number of nc-reads: " << reads_id.size() << endl;
    
    
    // get genome size
    size_t genome_size = get_genome_size(reads_range);
    
    // create a template to compare reads
    vector<bool> temp_array(genome_size*4+3, false);
    
    vector<int> cur_pu_var_count(4*(genome_size-1)+3+1, 0);
    vector<int> cur_pu_var_ref_count(4*(genome_size-1)+3+1, 0);
    vector<int> cur_pu_reads_count(genome_size, 0);
    
    // get topn nearest neighbors for each reads and find non-contained reads (no topn reads can cover all its range)
   
    int64_t n_nc_reads = 0;
    for (auto i : reads_id){
        ++n_nc_reads;
        if (n_nc_reads % 1000 == 0)
        cout << "processed " << n_nc_reads << " / " << reads_id.size() << endl;
        
        ConsensusSeq cur_cons;
        if (is_recode){
            cur_cons.seed = recode_data[i];
            cur_cons.cons_seq = recode_data[i];
        }else{
            cur_cons.seed = encode_data[i];
            cur_cons.cons_seq = encode_data[i];
        }
        cur_cons.start = reads_range[i].first;
        cur_cons.end = reads_range[i].second;
        
        bool is_valid = false;
        vector<int> prev_cons_seq = cur_cons.cons_seq;
        int b = 0;
        for (b = 0; b < max_iter; ++b){
            // calculate hamming distance between reads i and other reads and get topn nearest neighbors
            priority_queue<pair<int,double>, vector<pair<int,double> >, reads_compare_sim > topn_id;
            for (auto j = 0; j < recode_data.size(); ++j){
                if (i == j) continue;
                
                if (cur_cons.start >= reads_range[j].second || cur_cons.end <= reads_range[j].first)
                continue;
                
                if (cur_cons.start < reads_range[j].first && is_hang)
                continue;
                
                int min_overlap = int(0.75*(cur_cons.end - cur_cons.start));
                double cur_dist;
               
                ReadRange cur_range(cur_cons.start, cur_cons.end);
                if (is_recode){
                    cur_dist = sim_jaccard(cur_cons.cons_seq, recode_data[j], cur_range, reads_range[j], temp_array, true, min_overlap, true);
                }else{
                    cur_dist = sim_jaccard(cur_cons.cons_seq, encode_data[j], cur_range, reads_range[j], temp_array, true, min_overlap, true);
                }
                
                if (cur_dist <= min_jaccard) continue;
                
                topn_id.push(pair<int,double>(j,cur_dist));
            }
            
            // get max_nn neighbors
            vector<int> cur_neighbors;
            vector<double> cur_neighbors_dist;
            for (auto j = 0; j < max_nn; ++j){
                if (topn_id.empty()) break;
                
                int cur_id = topn_id.top().first;
                double cur_dist = topn_id.top().second;
                cur_neighbors.push_back(cur_id);
                cur_neighbors_dist.push_back(cur_dist);
        
                topn_id.pop();
            }
            
            if (cur_neighbors.size() < topn || cur_neighbors.size() < min_cvg)
            break;
            
            // pileup all neighbors and pop from most distant neighbor until all loci are homogeneous
            // pileup topn neighbors
            
            unordered_set<int64_t> mod_idx_var;
            unordered_set<int64_t> mod_idx_var_ref;
            for (auto j : cur_neighbors){
                pileup_var_online_count(cur_pu_var_count, recode_data[j], mod_idx_var);
                pileup_var_online_count(cur_pu_var_ref_count, recode_ref_data[j], mod_idx_var_ref);
            }
            
            // check pileup of all the neighbors
            bool is_homo = this->check_pileup_recode(cur_pu_var_count, cur_pu_var_ref_count, cur_cons.start, cur_cons.end, vector<int>(), min_cvg, min_prop, max_prop);
            
            // if not all loci are homogeneous, pop neighbors from the most distant one until all loci are homogeneous or number of neighbors <= topn
            int t = (int)cur_neighbors.size();
            if (!is_homo){
                for (auto j = 0; j < cur_neighbors.size(); ++j){
                    if (t <= topn) break;
                    t = (int)cur_neighbors.size() - 1 - j;
                    pileup_var_online_count_pop(cur_pu_var_count, recode_data[cur_neighbors[t]]);
                    pileup_var_online_count_pop(cur_pu_var_ref_count, recode_ref_data[cur_neighbors[t]]);
                    is_homo = this->check_pileup_recode(cur_pu_var_count, cur_pu_var_ref_count, cur_cons.start, cur_cons.end, vector<int>(), min_cvg, min_prop, max_prop);
                    if (is_homo) break;
                }
            }
            cur_cons.neighbors_id = vector<int> (cur_neighbors.begin(),  cur_neighbors.begin() + t);
            
            if (is_homo){
                this->get_consensus_recode(cur_cons, cur_pu_var_count, cur_pu_var_ref_count, cur_cons.start, cur_cons.end, min_cvg);
                is_valid = true;
                
                if (std::equal(cur_cons.cons_seq.begin(), cur_cons.cons_seq.end(), prev_cons_seq.begin())){
                    // clean cur_pu_var_count and cur_pu_reads_count
                    for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
                        cur_pu_var_count[*it] = 0;
                    for (auto it = mod_idx_var_ref.begin(); it != mod_idx_var_ref.end(); ++it)
                        cur_pu_var_ref_count[*it] = 0;
                    
                    break;
                }else{
                    prev_cons_seq = cur_cons.cons_seq;
                }
            }else{
                is_valid = false;
                
                // clean cur_pu_var_count and cur_pu_reads_count
                for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
                    cur_pu_var_count[*it] = 0;
                for (auto it = mod_idx_var_ref.begin(); it != mod_idx_var_ref.end(); ++it)
                    cur_pu_var_ref_count[*it] = 0;
                
                break;
            }
            
            // clean cur_pu_var_count and cur_pu_reads_count
            for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
                cur_pu_var_count[*it] = 0;
            for (auto it = mod_idx_var_ref.begin(); it != mod_idx_var_ref.end(); ++it)
                cur_pu_var_ref_count[*it] = 0;
        }
        if (is_valid){
            this->thread_locker.lock();
            rl_ann_clust.push_back(cur_cons);
            this->thread_locker.unlock();
        }
    }
}


void Assembler::print_rl_ann_clust(string outfile, bool is_metric, vector<int64_t> idx)
{
    if (idx.size() == 0){
        idx.resize(this->rl_ann_clust.size());
        for (auto i = 0; i < this->rl_ann_clust.size(); ++i)
            idx[i] = i;
    }
    
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    
    //for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
    for (auto &i : idx){
        if (rl_ann_clust[i].cons_seq.size() == 0)
            fs_outfile << -1 << '\t';
        else
            fs_outfile << this->rl_ann_clust[i].cons_seq << '\t';
        
        fs_outfile << this->rl_ann_clust[i].start << '\t' << this->rl_ann_clust[i].end << '\t';
        fs_outfile << this->rl_ann_clust[i].contig_count << '\t' << this->rl_ann_clust[i].contig_cvg << '\t';
        fs_outfile << this->rl_ann_clust[i].log_bf_null << '\t' << this->rl_ann_clust[i].log_bf_ind << '\t';
        fs_outfile << this->rl_ann_clust[i].rr_null << '\t' << this->rl_ann_clust[i].rr_ind;
        if (is_metric){
            fs_outfile << '\t';
            
            if (this->rl_ann_clust[i].seed.size() == 0)
                fs_outfile << -1 << '\t';
            else
                fs_outfile << this->rl_ann_clust[i].seed << '\t';
            
            if (rl_ann_clust[i].neighbors_id.size() == 0)
                fs_outfile << -1 << '\t';
            else
                fs_outfile << rl_ann_clust[i].neighbors_id << '\t';
            
            if (rl_ann_clust[i].tested_loci.size() == 0)
                fs_outfile << -1 << '\t';
            else
                fs_outfile << rl_ann_clust[i].tested_loci << '\t';
            
            if (rl_ann_clust[i].nn_reads_id.size() == 0)
                fs_outfile << -1;
            else
                fs_outfile << rl_ann_clust[i].nn_reads_id;
        }
        
        fs_outfile << endl;
    }
    
    
    fs_outfile.close();
        
}

void Assembler::print_nc_reads_id(string outfile)
{
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for (int64_t i = 0; i < this->nc_reads_id.size(); ++i)
        fs_outfile << nc_reads_id[i] << endl;
    fs_outfile.close();
}

void Assembler::read_ann_results(string ann_file,  bool is_ref)
{
    this->rl_ann_clust.clear();
    ifstream fs_ann_file;
    open_infile(fs_ann_file, ann_file);
    int64_t n_lines = 0;
    while(true){
        string buf;
        getline(fs_ann_file, buf);
        if(fs_ann_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size() != 12 && buf_vec.size() != 13){
            cerr << "buf_vec.size() = " << buf_vec.size() << endl;
            throw runtime_error("Line " + to_string(n_lines) + ": Assembler::read_ann_results, buf_vect.size() != 12 && buf_vec.size() != 13");
        }
        ConsensusSeq cur_cons;
        vector<int> cur_tmp_0 = split_int(buf_vec[0], ',');
        if (cur_tmp_0.size() >= 1){
            if (cur_tmp_0[0] != -1){
                cur_cons.cons_seq = cur_tmp_0;
            }else{
                if (is_ref){
                    cur_cons.cons_seq = vector<int>();
                }
            }
        }
        
        cur_cons.start = stoi(buf_vec[1]);
        cur_cons.end = stoi(buf_vec[2]);
        cur_cons.contig_count = stod(buf_vec[3]);
        cur_cons.contig_cvg = stod(buf_vec[4]);
        cur_cons.log_bf_null = stod(buf_vec[5]);
        cur_cons.log_bf_ind = stod(buf_vec[6]);
        cur_cons.rr_null = stod(buf_vec[7]);
        cur_cons.rr_ind = stod(buf_vec[8]);
        
        vector<int> cur_tmp_9 = split_int(buf_vec[9], ',');
        if (cur_tmp_9.size() >= 1){
            if (cur_tmp_9[0] != -1){
                cur_cons.seed = cur_tmp_9;
            }
        }

        vector<int> cur_tmp_10 = split_int(buf_vec[10], ',');
        if (cur_tmp_10.size() >= 1){
            if (cur_tmp_10[0] != -1){
                cur_cons.neighbors_id = cur_tmp_10;
            }
        }
        
        vector<int> cur_tmp_11 = split_int(buf_vec[11], ',');
        if (cur_tmp_11.size() >= 1){
            if (cur_tmp_11[0] != -1){
                cur_cons.tested_loci = cur_tmp_11;
            }
        }
                
        if (buf_vec.size() == 13){
            vector<int64_t> cur_tmp = split_int64_t(buf_vec[12], ',');
            if (cur_tmp.size() >= 1){
                if (cur_tmp[0] != -1){
                    cur_cons.nn_reads_id = cur_tmp;
                }
            }
            
        }
        
        this->rl_ann_clust.push_back(cur_cons);
        
        ++n_lines;
    }
    
    fs_ann_file.close();
    
}

void Assembler::filter_ann(string ann_file, double min_log_bf, double max_loci, double min_rr)
{
    this->read_ann_results(ann_file);
       
    vector<int64_t> idx_ft;
    
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        if ( ((rl_ann_clust[i].log_bf_null >= min_log_bf || rl_ann_clust[i].log_bf_ind >= min_log_bf) &&
            (rl_ann_clust[i].rr_null >= min_rr || rl_ann_clust[i].rr_ind >= min_rr)) ||
            rl_ann_clust[i].cons_seq.size() >= max_loci){
            idx_ft.push_back(i);
        }
    }
    
    this->print_rl_ann_clust(ann_file + ".ft", true, idx_ft);
    
    
}

void Assembler::ann_to_graph(Graph &gp, string ann_file, double min_prop, double min_len_prop, double min_jaccard, bool is_or, int max_ndiff)
{
    this->read_ann_results(ann_file);
    if (this->rl_ann_clust.size() == 0)
        return;
    
    size_t gsize = 0;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        add_vertex(gp);
        if (this->rl_ann_clust[i].end > gsize)
            gsize = this->rl_ann_clust[i].end;
    }
    vector<bool> temp_array(gsize*4 + 3, false);
    
    // get maximal encoded cons_seq and tested_loci
    int64_t temp_size = 0;
    int64_t temp_tested_size = 0;
    for (int64_t i = 0; i < this->rl_ann_clust.size(); ++i){
        for (int64_t j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j)
            if (rl_ann_clust[i].cons_seq[j] > temp_size)
                temp_size = rl_ann_clust[i].cons_seq[j];
        
        for (int64_t j = 0; j < this->rl_ann_clust[i].tested_loci.size(); ++j)
            if (rl_ann_clust[i].tested_loci[j] > temp_tested_size)
                temp_tested_size = rl_ann_clust[i].tested_loci[j];
        
    }
    if (temp_size == 0)
        throw runtime_error("Assembler::ann_to_graph, temp_size == 0");
    if (temp_tested_size == 0)
        throw runtime_error("Assembler::ann_to_graph, temp_tested_size == 0");
    
    ++temp_size;
    ++temp_tested_size;
    
    // generate a template vector, note that temp_vec is as least 4 times larger than temp_tested_vec
    vector<bool> temp_vec_i(temp_size, false);
    vector<bool> temp_vec_j(temp_size, false);
    vector<bool> temp_tested_vec_i(temp_tested_size, false);
    vector<bool> temp_tested_vec_j(temp_tested_size, false);
    
    // pairwise compare contigs
    for (int64_t i = 0; i < rl_ann_clust.size(); ++i){
        // fill in template by the ith cons_seq and tested_loci
        for (int64_t k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k)
            temp_vec_i[rl_ann_clust[i].cons_seq[k]] = true;
        
        for (int64_t k = 0; k < rl_ann_clust[i].tested_loci.size(); ++k)
            temp_tested_vec_i[rl_ann_clust[i].tested_loci[k]] = true;
        
        // compare contigs
        //int start_code = 4*rl_ann_clust[i].start;
        //int end_code = 4*rl_ann_clust[i].end +3;
        
        for (int64_t j = 0; j < this->rl_ann_clust.size(); ++j){
            if (j == i) continue;

            if (rl_ann_clust[j].start < rl_ann_clust[i].start || rl_ann_clust[j].start > rl_ann_clust[i].end || rl_ann_clust[j].end < rl_ann_clust[i].end)
                continue;
            
            if (rl_ann_clust[j].start == rl_ann_clust[i].start && rl_ann_clust[j].end == rl_ann_clust[i].end && j > i)
                continue;
            
            int overlap_len = rl_ann_clust[i].end - rl_ann_clust[j].start + 1;
            int contig_len_i = rl_ann_clust[i].end - rl_ann_clust[i].start + 1;
            int contig_len_j = rl_ann_clust[j].end - rl_ann_clust[j].start + 1;
            
            //bool is_len_overlap_ok =
            /*if (is_or){
                if (overlap_len < min_len_prop*(rl_ann_clust[i].end - rl_ann_clust[i].start + 1) &&
                    overlap_len < min_len_prop*(rl_ann_clust[j].end - rl_ann_clust[j].start + 1))
                    continue;
            }else{
                if (overlap_len < min_len_prop*(rl_ann_clust[i].end - rl_ann_clust[i].start + 1) ||
                    overlap_len < min_len_prop*(rl_ann_clust[j].end - rl_ann_clust[j].start + 1))
                    continue;
            }*/
            
            double n_cons_seq_i = rl_ann_clust[i].cons_seq.size();
            //double n_cons_seq_j = 0;
            double n_cons_seq_j = rl_ann_clust[j].cons_seq.size();
            double n_overlap = 0;
            
            // fill in template by the jth cons_seq and tested_loci
            bool is_diff = false;
            int64_t n_diff = 0;
            for (int64_t k = 0; k < rl_ann_clust[j].cons_seq.size(); ++k){
                temp_vec_j[rl_ann_clust[j].cons_seq[k]] = true;
                if (!temp_vec_i[rl_ann_clust[j].cons_seq[k]] && temp_tested_vec_i[rl_ann_clust[j].cons_seq[k] / 4]){
                    is_diff = true;
                    ++n_diff;
                }
                
                if (temp_vec_i[rl_ann_clust[j].cons_seq[k]])
                    ++n_overlap;
                
                //if (rl_ann_clust[j].cons_seq[k] >= start_code && rl_ann_clust[j].cons_seq[k] <= end_code)
                //    ++n_cons_seq_j;
            }
            
            for (int64_t k = 0; k < rl_ann_clust[j].tested_loci.size(); ++k){
                temp_tested_vec_j[rl_ann_clust[j].tested_loci[k]] = true;
            }
            
            // scan the ith contig again
            for (int64_t k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k){
                if (!temp_vec_j[rl_ann_clust[i].cons_seq[k]] && temp_tested_vec_j[rl_ann_clust[i].cons_seq[k] / 4]){
                    is_diff = true;
                    ++n_diff;
                }
            }
            
            double cur_jaccard = sim_jaccard(rl_ann_clust[i].cons_seq, rl_ann_clust[j].cons_seq,
                                             ReadRange(rl_ann_clust[i].start, rl_ann_clust[i].end),
                                             ReadRange(rl_ann_clust[j].start, rl_ann_clust[j].end),
                                             temp_array, false, 0, false);
            //cout << i << ',' << j << " : " << cur_jaccard << endl;
            
            if (is_or){
                if ((!is_diff || (cur_jaccard >= min_jaccard && n_diff <= max_ndiff)) && (n_overlap > min_prop*n_cons_seq_i || n_overlap > min_prop*n_cons_seq_j || overlap_len > min_len_prop*contig_len_i || overlap_len > min_len_prop*contig_len_j ))
                //if ((!is_diff || (cur_jaccard >= min_jaccard || n_diff <= 1)) && (n_overlap > min_prop*n_cons_seq_i || n_overlap > min_prop*n_cons_seq_j || overlap_len > min_len_prop*contig_len_i || overlap_len > min_len_prop*contig_len_j ))
                    boost::add_edge(i, j, gp);
            }else{
                if ((!is_diff || (cur_jaccard >= min_jaccard && n_diff <= max_ndiff)) && (n_overlap > min_prop*n_cons_seq_i && n_overlap > min_prop*n_cons_seq_j && overlap_len > min_len_prop*contig_len_i && overlap_len > min_len_prop*contig_len_j ))
                    boost::add_edge(i, j, gp);
            }
            
            // clear template of the jth contig
            for (int64_t k = 0; k < rl_ann_clust[j].cons_seq.size(); ++k){
                temp_vec_j[rl_ann_clust[j].cons_seq[k]] = false;
            }
            
            for (int64_t k = 0; k < rl_ann_clust[j].tested_loci.size(); ++k){
                temp_tested_vec_j[rl_ann_clust[j].tested_loci[k]] = false;
            }
            
        }
        
        // clear template
        for (int64_t k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k)
            temp_vec_i[rl_ann_clust[i].cons_seq[k]] = false;
        
        for (int64_t k = 0; k < rl_ann_clust[i].tested_loci.size(); ++k)
            temp_tested_vec_i[rl_ann_clust[i].tested_loci[k]] = false;
        
    }
    
}

void Assembler::test_contigs_refonly(const vector<vector<int> > &recode_data, const vector<vector<int> > &recode_ref_data, const vector<ReadRange> &reads_range)
{
    if (recode_data.size() != recode_ref_data.size() || recode_data.size()!= reads_range.size())
        throw runtime_error("incompatible recode_data, recode_ref_data or reads_range");
    
    // check homo_blocks
    int64_t max_range = -1;
    for (auto i = 0; i < reads_range.size(); ++i)
        if (reads_range[i].second > max_range)
            max_range = reads_range[i].second;
    if (max_range + 1 > this->homo_blocks.size())
        throw runtime_error("Assembler::test_contigs, max_range + 1 > this->homo_blocks.size().");
    
    vector<vector<int> > pu_recode_full = pileup_var(recode_data);
    vector<vector<int> > pu_recode_ref_full = pileup_var(recode_ref_data);
    
    
    // test single locus contig
    ConsensusSeq ann_ref = this->rl_ann_clust.back();
    if (ann_ref.cons_seq.size() != 0)
        throw runtime_error(" Assembler::test_contigs() : ann_ref.cons_seq.size() != 0");
    
    cout << "test single locus contigs" << endl;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        // get reads on contig i and ref cons
        unordered_set<int64_t> nn_reads_ids(this->rl_ann_clust[i].nn_reads_id.begin(), this->rl_ann_clust[i].nn_reads_id.end());
        nn_reads_ids.insert(ann_ref.nn_reads_id.begin(), ann_ref.nn_reads_id.end());
        
        unordered_map<int64_t, vector<int> > pu_recode;
        unordered_map<int64_t, vector<int> > pu_recode_ref;
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j){
            int64_t cur_code = this->rl_ann_clust[i].cons_seq[j];
            int64_t cur_locus = cur_code / 4;
            
            for (auto k = 0; k < 4; ++k){
                for (auto cur_readid : pu_recode_full[4*cur_locus + k]){
                    if (nn_reads_ids.find(cur_readid) != nn_reads_ids.end())
                        pu_recode[4*cur_locus + k].push_back(cur_readid);
                }
                
                for (auto cur_readid : pu_recode_ref_full[4*cur_locus + k]){
                    if (nn_reads_ids.find(cur_readid) != nn_reads_ids.end())
                        pu_recode_ref[4*cur_locus + k].push_back(cur_readid);
                }
                //pu_recode[4*cur_locus + k] = pu_recode_full[4*cur_locus + k];
                //pu_recode_ref[4*cur_locus + k] =  pu_recode_ref_full[4*cur_locus + k];
            }
        }
        
        // test contigs
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j){
            int cur_count = 0;
            int cur_cvg = 0;
            
            int64_t cur_code = this->rl_ann_clust[i].cons_seq[j];
            int64_t cur_locus = cur_code / 4;
            
            // calculate marginal probability
            cur_count = (int)pu_recode[cur_code].size();
            cur_cvg = (int)pu_recode[4*cur_locus].size() + (int)pu_recode[4*cur_locus + 1].size() + (int)pu_recode[4*cur_locus + 2].size() + (int)pu_recode[4*cur_locus + 3].size();
            cur_cvg += pu_recode_ref[4*cur_locus].size() + pu_recode_ref[4*cur_locus + 1].size() + pu_recode_ref[4*cur_locus + 2].size() + pu_recode_ref[4*cur_locus + 3].size();
            
            if (cur_cvg > 0 ){
                if (this->alpha == -1 || this->beta == -1)
                    throw runtime_error("this->alpha == -1 || this->beta == -1");
                
                double cur_log_bf_null = binom_log_bf(cur_count, cur_cvg, this->alpha, this->beta);
                if (cur_log_bf_null < this->rl_ann_clust[i].log_bf_null || j == 0)
                    this->rl_ann_clust[i].log_bf_null = cur_log_bf_null;
                
                double cur_rr = cur_count*(this->alpha + this->beta) / (cur_cvg*this->alpha);
                if (cur_rr < this->rl_ann_clust[i].rr_null || j ==0)
                    this->rl_ann_clust[i].rr_null = cur_rr;
            }
        }
    }
    
    // test independence
    cout << "test multiple loci contigs" << endl;
    
    vector<int64_t> temp_read_id(recode_data.size(), -1);
    vector<int64_t> temp_read_id_cvg(recode_data.size(), -1);
    
    vector<int64_t> temp_read_id_block(recode_data.size(), -1);
    vector<int64_t> temp_read_id_cvg_block(recode_data.size(), -1);
    
    int64_t counter = 0;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        //cout << i << endl;
        if (this->rl_ann_clust[i].cons_seq.size() < 2)
            continue;
        
        // get reads on contig i and ref cons
        unordered_set<int64_t> nn_reads_ids(this->rl_ann_clust[i].nn_reads_id.begin(), this->rl_ann_clust[i].nn_reads_id.end());
        nn_reads_ids.insert(ann_ref.nn_reads_id.begin(), ann_ref.nn_reads_id.end());
        
        unordered_map<int64_t, vector<int> > pu_recode;
        unordered_map<int64_t, vector<int> > pu_recode_ref;
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j){
            int64_t cur_code = this->rl_ann_clust[i].cons_seq[j];
            int64_t cur_locus = cur_code / 4;
            
            for (auto k = 0; k < 4; ++k){
                for (auto cur_readid : pu_recode_full[4*cur_locus + k]){
                    if (nn_reads_ids.find(cur_readid) != nn_reads_ids.end())
                        pu_recode[4*cur_locus + k].push_back(cur_readid);
                }
                
                for (auto cur_readid : pu_recode_ref_full[4*cur_locus + k]){
                    if (nn_reads_ids.find(cur_readid) != nn_reads_ids.end())
                        pu_recode_ref[4*cur_locus + k].push_back(cur_readid);
                }
                //pu_recode[4*cur_locus + k] = pu_recode_full[4*cur_locus + k];
                //pu_recode_ref[4*cur_locus + k] =  pu_recode_ref_full[4*cur_locus + k];
            }
        }
        
        // merge loci within MIN_LOCI_BLOCK_SIZE, and test if joint_prop > product(block_prop)
        vector<vector<int> > blocks;
        vector<double> block_prop, block_count, block_cvg;
        double cur_block_prop = -1, cur_block_count = 0, cur_block_cvg = 0;
        double joint_prop = -1, joint_count = 0, joint_cvg = 0;
        
        int64_t prev_locus = -MIN_LOCI_BLOCK_SIZE - 1;
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j){
            int64_t cur_code = this->rl_ann_clust[i].cons_seq[j];
            int64_t cur_locus = cur_code / 4;
            
            if (cur_locus < prev_locus)
                throw runtime_error("Assembler::test_contigs, cur_locus < prev_locus (locus in contigs is not sorted)");
            
            // check if enter a new block
            bool is_new_block = false;
            if (j == 0){
                is_new_block = true;
            }else{
                if (this->homo_blocks[cur_locus] - this->homo_blocks[prev_locus] > MIN_LOCI_BLOCK_SIZE)
                    is_new_block = true;
            }
            
            // record previous block info
            if (is_new_block){
                if (j > 0){
                    block_count.push_back(cur_block_count);
                    block_cvg.push_back(cur_block_cvg);
                    block_prop.push_back(cur_block_prop);
                }
                blocks.push_back(vector<int>());
            }
            blocks.back().push_back((int)cur_locus);
            
            // reset current count
            cur_block_prop = -1; cur_block_count = 0; cur_block_cvg = 0;
            joint_prop = -1; joint_count = 0; joint_cvg = 0;
            
            // fill template and count
            for (auto k = 0; k < pu_recode[cur_code].size(); ++k){
                // count joint probability
                if (temp_read_id[pu_recode[cur_code][k]] == counter - 1 || j == 0){
                    ++joint_count;
                    temp_read_id[pu_recode[cur_code][k]] = counter;
                }
                
                // count block probability
                if (temp_read_id_block[pu_recode[cur_code][k]] == counter - 1 || is_new_block){
                    ++cur_block_count;
                    temp_read_id_block[pu_recode[cur_code][k]] = counter;
                }
                
            }
            
            for (auto s = 0; s <= 3; ++s){
                for (auto k = 0; k < pu_recode[4*cur_locus + s].size(); ++k){
                    // count joint probability
                    if (temp_read_id_cvg[pu_recode[4*cur_locus + s][k]] == counter - 1 || j == 0){
                        ++joint_cvg;
                        temp_read_id_cvg[pu_recode[4*cur_locus + s][k]] = counter;
                    }
                    
                    // count block probability
                    if (temp_read_id_cvg_block[pu_recode[4*cur_locus + s][k]] == counter - 1 || is_new_block){
                        ++cur_block_cvg;
                        temp_read_id_cvg_block[pu_recode[4*cur_locus + s][k]] = counter;
                    }
                }
                
                for (auto k = 0; k < pu_recode_ref[4*cur_locus + s].size(); ++k){
                    // count joint probability
                    if (temp_read_id_cvg[pu_recode_ref[4*cur_locus + s][k]] == counter - 1 || j == 0){
                        ++joint_cvg;
                        temp_read_id_cvg[pu_recode_ref[4*cur_locus + s][k]] = counter;
                    }
                    
                    // count block probability
                    if (temp_read_id_cvg_block[pu_recode_ref[4*cur_locus + s][k]] == counter - 1 || is_new_block){
                        ++cur_block_cvg;
                        temp_read_id_cvg_block[pu_recode_ref[4*cur_locus + s][k]] = counter;
                    }
                }
                
            }
            
            ++counter;
            if (counter >= numeric_limits<int64_t>::max() - 1)
                throw runtime_error("Assembler::test_contigs, counter >= numeric_limits<int64_t>::max() - 1");
            
            // calculate block probability
            if (cur_block_cvg > 0)
                cur_block_prop = cur_block_count / cur_block_cvg;
            
            // calculater current joint probablity
            if (joint_cvg > 0 && j == this->rl_ann_clust[i].cons_seq.size() - 1)
                joint_prop = joint_count / joint_cvg;
            
            // record block info
            if (j == this->rl_ann_clust[i].cons_seq.size() - 1){
                block_count.push_back(cur_block_count);
                block_cvg.push_back(cur_block_cvg);
                block_prop.push_back(cur_block_prop);
            }
            
            prev_locus = cur_locus;
        }
        
        // don't test if joint count is 0 (two many loci so too few reads having all of them)
        if (joint_count == 0)
            continue;
        
        // only test independence if there are two or more
        if (block_prop.size() >= 2){
            double exp_prop = prod(block_prop);
            if (exp_prop > 0){
                this->rl_ann_clust[i].log_bf_ind = binom_log_bf(joint_count, joint_cvg, exp_prop);
                this->rl_ann_clust[i].rr_ind = joint_count / (joint_cvg*exp_prop);
                
            }
        }
        
    }
}

void Assembler::test_contigs(const vector<vector<int> > &recode_data, const vector<vector<int> > &recode_ref_data, const vector<ReadRange> &reads_range)
{
    if (recode_data.size() != recode_ref_data.size() || recode_data.size()!= reads_range.size())
        throw runtime_error("incompatible recode_data, recode_ref_data or reads_range");

    // check homo_blocks
    int64_t max_range = -1;
    for (auto i = 0; i < reads_range.size(); ++i)
        if (reads_range[i].second > max_range)
            max_range = reads_range[i].second;
    if (max_range + 1 > this->homo_blocks.size())
        throw runtime_error("Assembler::test_contigs, max_range + 1 > this->homo_blocks.size().");
        
    vector<vector<int> > pu_recode = pileup_var(recode_data);
    vector<vector<int> > pu_recode_ref = pileup_var(recode_ref_data);
    
    // test single locus contig
    cout << "test single locus contigs" << endl;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j){
            int cur_count = 0;
            int cur_cvg = 0;
            
            int64_t cur_code = this->rl_ann_clust[i].cons_seq[j];
            int64_t cur_locus = cur_code / 4;
            
            // calculate marginal probability
            cur_count = (int)pu_recode[cur_code].size();
            cur_cvg = (int)pu_recode[4*cur_locus].size() + (int)pu_recode[4*cur_locus + 1].size() + (int)pu_recode[4*cur_locus + 2].size() + (int)pu_recode[4*cur_locus + 3].size();
            cur_cvg += pu_recode_ref[4*cur_locus].size() + pu_recode_ref[4*cur_locus + 1].size() + pu_recode_ref[4*cur_locus + 2].size() + pu_recode_ref[4*cur_locus + 3].size();
            
            if (cur_cvg > 0 ){
                if (this->alpha == -1 || this->beta == -1)
                    throw runtime_error("this->alpha == -1 || this->beta == -1");
                
                double cur_log_bf_null = binom_log_bf(cur_count, cur_cvg, this->alpha, this->beta);
                if (cur_log_bf_null < this->rl_ann_clust[i].log_bf_null || j == 0)
                    this->rl_ann_clust[i].log_bf_null = cur_log_bf_null;
                
                double cur_rr = cur_count*(this->alpha + this->beta) / (cur_cvg*this->alpha);
                if (cur_rr < this->rl_ann_clust[i].rr_null || j ==0)
                    this->rl_ann_clust[i].rr_null = cur_rr;
            }
        }
    }
    
    // test independence
    cout << "test multiple loci contigs" << endl;
    
    vector<int64_t> temp_read_id(recode_data.size(), -1);
    vector<int64_t> temp_read_id_cvg(recode_data.size(), -1);
    
    vector<int64_t> temp_read_id_block(recode_data.size(), -1);
    vector<int64_t> temp_read_id_cvg_block(recode_data.size(), -1);
    
    int64_t counter = 0;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        //cout << i << endl;
        if (this->rl_ann_clust[i].cons_seq.size() < 2)
            continue;
        
        // merge loci within MIN_LOCI_BLOCK_SIZE, and test if joint_prop > product(block_prop)
        vector<vector<int> > blocks;
        vector<double> block_prop, block_count, block_cvg;
        double cur_block_prop = -1, cur_block_count = 0, cur_block_cvg = 0;
        double joint_prop = -1, joint_count = 0, joint_cvg = 0;
        
        int64_t prev_locus = -MIN_LOCI_BLOCK_SIZE - 1;
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j){
            int64_t cur_code = this->rl_ann_clust[i].cons_seq[j];
            int64_t cur_locus = cur_code / 4;
    
            if (cur_locus < prev_locus)
                throw runtime_error("Assembler::test_contigs, cur_locus < prev_locus (locus in contigs is not sorted)");
            
            // check if enter a new block
            bool is_new_block = false;
            if (j == 0){
                is_new_block = true;
            }else{
                if (this->homo_blocks[cur_locus] - this->homo_blocks[prev_locus] > MIN_LOCI_BLOCK_SIZE)
                    is_new_block = true;
            }
            
            // record previous block info
            if (is_new_block){
                if (j > 0){
                    block_count.push_back(cur_block_count);
                    block_cvg.push_back(cur_block_cvg);
                    block_prop.push_back(cur_block_prop);
                }
                blocks.push_back(vector<int>());
            }
            blocks.back().push_back((int)cur_locus);
            
            // reset current count
            cur_block_prop = -1; cur_block_count = 0; cur_block_cvg = 0;
            joint_prop = -1; joint_count = 0; joint_cvg = 0;
            
            // fill template and count
            for (auto k = 0; k < pu_recode[cur_code].size(); ++k){
                // count joint probability
                if (temp_read_id[pu_recode[cur_code][k]] == counter - 1 || j == 0){
                    ++joint_count;
                    temp_read_id[pu_recode[cur_code][k]] = counter;
                }
                
                // count block probability
                if (temp_read_id_block[pu_recode[cur_code][k]] == counter - 1 || is_new_block){
                    ++cur_block_count;
                    temp_read_id_block[pu_recode[cur_code][k]] = counter;
                }
                
            }
            
            for (auto s = 0; s <= 3; ++s){
                for (auto k = 0; k < pu_recode[4*cur_locus + s].size(); ++k){
                    // count joint probability
                    if (temp_read_id_cvg[pu_recode[4*cur_locus + s][k]] == counter - 1 || j == 0){
                        ++joint_cvg;
                        temp_read_id_cvg[pu_recode[4*cur_locus + s][k]] = counter;
                    }
                    
                    // count block probability
                    if (temp_read_id_cvg_block[pu_recode[4*cur_locus + s][k]] == counter - 1 || is_new_block){
                        ++cur_block_cvg;
                        temp_read_id_cvg_block[pu_recode[4*cur_locus + s][k]] = counter;
                    }
                }
                
                for (auto k = 0; k < pu_recode_ref[4*cur_locus + s].size(); ++k){
                    // count joint probability
                    if (temp_read_id_cvg[pu_recode_ref[4*cur_locus + s][k]] == counter - 1 || j == 0){
                        ++joint_cvg;
                        temp_read_id_cvg[pu_recode_ref[4*cur_locus + s][k]] = counter;
                    }
                    
                    // count block probability
                    if (temp_read_id_cvg_block[pu_recode_ref[4*cur_locus + s][k]] == counter - 1 || is_new_block){
                        ++cur_block_cvg;
                        temp_read_id_cvg_block[pu_recode_ref[4*cur_locus + s][k]] = counter;
                    }
                }
                
            }
            
            ++counter;
            if (counter >= numeric_limits<int64_t>::max() - 1)
                throw runtime_error("Assembler::test_contigs, counter >= numeric_limits<int64_t>::max() - 1");
            
            // calculate block probability
            if (cur_block_cvg > 0)
                cur_block_prop = cur_block_count / cur_block_cvg;
            
            // calculater current joint probablity
            if (joint_cvg > 0 && j == this->rl_ann_clust[i].cons_seq.size() - 1)
                joint_prop = joint_count / joint_cvg;
            
            // record block info
            if (j == this->rl_ann_clust[i].cons_seq.size() - 1){
                block_count.push_back(cur_block_count);
                block_cvg.push_back(cur_block_cvg);
                block_prop.push_back(cur_block_prop);
            }

            prev_locus = cur_locus;
        }
        
        // don't test if joint count is 0 (two many loci so too few reads having all of them)
        if (joint_count == 0)
            continue;
        
        // only test independence if there are two or more
        if (block_prop.size() >= 2){
            double exp_prop = prod(block_prop);
            if (exp_prop > 0){
                this->rl_ann_clust[i].log_bf_ind = binom_log_bf(joint_count, joint_cvg, exp_prop);
                this->rl_ann_clust[i].rr_ind = joint_count / (joint_cvg*exp_prop);
                
            }
        }
        
    }
    
}

void Assembler::test_contigs_pairwise(string ann_file, string recode_file, string out_file, double min_log_bf, int max_loci, int min_cvg, double min_rr)
{
    cout << "load ann" << endl;
    this->read_ann_results(ann_file);
    
    cout << "load recode and recode_ref" << endl;
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_file + ".ref");
    
    if (recode_data.size() != recode_ref_data.size())
        throw runtime_error("recode_data.size() != recode_ref_data.size()");
    
    cout << "pileup" << endl;
    vector<vector<int> > pu_recode = pileup_var(recode_data);
    vector<vector<int> > pu_recode_ref = pileup_var(recode_ref_data);
    
    if (pu_recode.size()/4 + 1 > this->homo_blocks.size() || pu_recode_ref.size()/4 + 1 > this->homo_blocks.size())
        throw runtime_error("Assembler::test_contigs_pairwise, pu_recode.size() + 1 > this->homo_blocks.size() || pu_recode_ref.size() + 1 > this->homo_blocks.size().");
    
    // get maximal contig locus
    size_t contig_size = 0;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i)
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j)
            if (this->rl_ann_clust[i].cons_seq[j] + 1 > contig_size)
                contig_size = this->rl_ann_clust[i].cons_seq[j] + 1;
    if (contig_size == 0)
        return;
    
    // merge close contigs
    vector<bool> temp_reads(recode_data.size(), false);
    vector<int> temp_var(contig_size, 0);
    vector<int64_t> temp_count(recode_data.size(), -1);
    vector<int64_t> temp_count_cvg(recode_data.size(), -1);
    int64_t counter = -1;
    vector<int64_t> idx;
    
    cout << "test contigs:" << endl;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        if ((i+1) % 100 == 0) cout << "tested " << i+1 << " contigs" << endl;
        
        if (rl_ann_clust[i].cons_seq.size() == 0)
            continue;
        
        // fill temp_var
        for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k)
            temp_var[this->rl_ann_clust[i].cons_seq[k]] = 1;
        
        // get close neighbor contigs
        vector<int64_t> neighbor_ids;
        bool is_noise = false;
        
        for (auto j = 0; j < this->rl_ann_clust.size(); ++j){
            if (j == i) continue;
            
            vector<int> diff_var;
            vector<bool> diff_var_type; // true for putative false positives (i has j doesn't) and false for putative false negatives (j has i doesn't)
            
            // check overlap of two contigs
            int64_t overlap_start = this->rl_ann_clust[i].start >= this->rl_ann_clust[j].start ?
                                    this->rl_ann_clust[i].start : this->rl_ann_clust[j].start;
            
            int64_t overlap_end = this->rl_ann_clust[i].end <= this->rl_ann_clust[j].end ?
                                  this->rl_ann_clust[i].end : this->rl_ann_clust[j].end;
            
            int64_t overlap_len = overlap_end - overlap_start + 1;
            
            if(overlap_len <= 0) continue;
            
            vector<int> cons_seq_i;
            vector<int> cons_seq_j;
            
            for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k){
                if (this->rl_ann_clust[i].cons_seq[k] >= 4*overlap_start &&
                    this->rl_ann_clust[i].cons_seq[k] <= 4*overlap_end + 3){
                    cons_seq_i.push_back(this->rl_ann_clust[i].cons_seq[k]);
                }
            }
            
            for (auto k = 0; k < this->rl_ann_clust[j].cons_seq.size(); ++k){
                if (this->rl_ann_clust[j].cons_seq[k] >= 4*overlap_start &&
                    this->rl_ann_clust[j].cons_seq[k] <= 4*overlap_end + 3){
                    cons_seq_j.push_back(this->rl_ann_clust[j].cons_seq[k]);
                }
            }
            
            if (cons_seq_i.size() < 0.75 * this->rl_ann_clust[i].cons_seq.size())
                continue;
            
            // check similarity
            int64_t n_common = 0;
            int64_t n_union = cons_seq_i.size();
            
            for (auto k = 0; k < cons_seq_j.size(); ++k){
                if (temp_var[cons_seq_j[k]] == 1){
                    ++n_common;
                    temp_var[cons_seq_j[k]] = 2;
                }else{
                    ++n_union;
                }
            }
            
            int64_t n_diff = n_union - n_common;
            if ( n_diff > max_loci || this->rl_ann_clust[i].contig_cvg > this->rl_ann_clust[j].contig_cvg ||
                ( n_common < int(0.5*(cons_seq_i.size()))) ){
                
                for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k)
                    temp_var[this->rl_ann_clust[i].cons_seq[k]] = 1;
                continue;
            }
            
            // record different variants
            neighbor_ids.push_back(j);
            
            unordered_set<int> cur_tested_loci_i(this->rl_ann_clust[i].tested_loci.begin(), this->rl_ann_clust[i].tested_loci.end());
            unordered_set<int> cur_tested_loci_j(this->rl_ann_clust[j].tested_loci.begin(), this->rl_ann_clust[j].tested_loci.end());
            
            for (auto k = 0; k < cons_seq_i.size(); ++k){
                if (temp_var[cons_seq_i[k]] == 1 &&
                    cur_tested_loci_j.find(int(cons_seq_i[k] / 4)) != cur_tested_loci_j.end()){
                    diff_var.push_back(cons_seq_i[k]);
                    diff_var_type.push_back(true);
                }
            }
            
            for (auto k = 0; k < cons_seq_j.size(); ++k){
                if (temp_var[cons_seq_j[k]] == 0 &&
                    cur_tested_loci_i.find(int(cons_seq_j[k] / 4)) != cur_tested_loci_i.end()){
                    diff_var.push_back(cons_seq_j[k]);
                    diff_var_type.push_back(false);
                }
            }
            
            for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k)
                temp_var[this->rl_ann_clust[i].cons_seq[k]] = 1;
            
            
            // calculate joint frequency of the differences
            unordered_set<int64_t> nn_reads_ids(this->rl_ann_clust[i].nn_reads_id.begin(), this->rl_ann_clust[i].nn_reads_id.end());
            nn_reads_ids.insert(this->rl_ann_clust[j].nn_reads_id.begin(), this->rl_ann_clust[j].nn_reads_id.end());
            
            for (auto it = nn_reads_ids.begin(); it != nn_reads_ids.end(); ++it)
                temp_reads[*it] = true;
            
            double joint_count = 0;
            double joint_cvg = 0;
            for (auto s = 0; s < diff_var.size(); ++s){
                ++counter;
                if (counter >= numeric_limits<int64_t>::max() - 1)
                    throw runtime_error("Assembler::test_contigs, counter >= numeric_limits<int64_t>::max() - 1");
                
                int cur_locus = diff_var[s] / 4;
                joint_count = 0;
                joint_cvg = 0;
                
                // calculate joint count
                if (diff_var_type[s]){
                    for (auto k = 0; k < pu_recode[diff_var[s]].size(); ++k){
                        if ( temp_reads[pu_recode[diff_var[s]][k]] && (temp_count[pu_recode[diff_var[s]][k]] == counter - 1 || s == 0) ){
                            temp_count[pu_recode[diff_var[s]][k]] = counter;
                            ++joint_count;
                        }
                    }
                }else{
                    for (int t = 0; t <= 3; ++t){
                        for (auto k = 0; k < pu_recode_ref[4*cur_locus+t].size(); ++k){
                            if ( temp_reads[pu_recode_ref[4*cur_locus+t][k]] && (temp_count[pu_recode_ref[4*cur_locus+t][k]] == counter - 1 || s == 0) ){
                                temp_count[pu_recode_ref[4*cur_locus+t][k]] = counter;
                                ++joint_count;
                            }
                        }
                    }
                }
                
                // calculate joint coverage
                for (int t = 0; t <= 3; ++t){
                    for (auto k = 0; k < pu_recode[4*cur_locus+t].size(); ++k){
                        if ( temp_reads[ pu_recode[4*cur_locus+t][k] ] && ( temp_count_cvg[ pu_recode[4*cur_locus+t][k] ] == counter - 1 || s == 0) ){
                            temp_count_cvg[ pu_recode[4*cur_locus+t][k] ] = counter;
                            ++joint_cvg;
                        }
                    }
                    
                    for (auto k = 0; k < pu_recode_ref[4*cur_locus+t].size(); ++k){
                        if ( temp_reads[ pu_recode_ref[4*cur_locus+t][k] ] && ( temp_count_cvg[ pu_recode_ref[4*cur_locus+t][k] ] == counter - 1 || s == 0) ){
                            temp_count_cvg[ pu_recode_ref[4*cur_locus+t][k] ] = counter;
                            ++joint_cvg;
                        }
                    }
                }
            }
            
            for (auto it = nn_reads_ids.begin(); it != nn_reads_ids.end(); ++it)
                temp_reads[*it] = false;
            
            // calculate bayes factor
            sort(diff_var.begin(), diff_var.end());
            int n_blocks = 0;
            for (auto k = 0; k < diff_var.size(); ++k){
                if (k == 0){
                    ++n_blocks;
                }else{
                    if (this->homo_blocks[diff_var[k]/4] - this->homo_blocks[diff_var[k-1]/4] > MIN_LOCI_BLOCK_SIZE )
                        ++n_blocks;
                }
            }
            
            //if (n_blocks == 0) throw runtime_error("n_blocks == 0");
            if (n_blocks > 0){
                //cout << i << "," << j << endl;
                if (this->alpha == -1 || this->beta == -1)
                    throw runtime_error("this->alpha == -1 || this->beta == -1");
                
                double exp_prop = pow(this->alpha/(this->alpha + this->beta), n_blocks);
                if (exp_prop < EPS) throw runtime_error("exp_prop < EPS");
                
                if (joint_cvg >= min_cvg){
                    double cur_log_bf = binom_log_bf(joint_count, joint_cvg, exp_prop);
                    double cur_rr = joint_count / (joint_cvg * exp_prop);
                    if (cur_log_bf < min_log_bf || cur_rr < min_rr){
                        is_noise = true;
                        break;
                    }
                }
            }
        }
        
        // clean temp_var
        for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k)
            temp_var[this->rl_ann_clust[i].cons_seq[k]] = 0;
        
        if (!is_noise) idx.push_back(i);
    }
    
    // print results
    this->print_rl_ann_clust(out_file, true, idx);
}

void Assembler::test_contigs_pairwise_legacy(string ann_file, string recode_file, string out_file, double min_log_bf, int max_loci, int min_cvg)
{
    // load ann
    this->read_ann_results(ann_file);
    
    // load recode and recode_ref
    vector<vector<int> > recode_data;
    loadencodedata(recode_data, recode_file);
    
    vector<vector<int> > recode_ref_data;
    loadencodedata(recode_ref_data, recode_file + ".ref");
    
    if (recode_data.size() != recode_ref_data.size())
        throw runtime_error("recode_data.size() != recode_ref_data.size()");
        
    // pileup
    vector<vector<int> > pu_recode = pileup_var(recode_data);
    vector<vector<int> > pu_recode_ref = pileup_var(recode_ref_data);
    
    // get maximal contig locus
    size_t contig_size = 0;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i)
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j)
            if (this->rl_ann_clust[i].cons_seq[j] + 1 > contig_size)
                contig_size = this->rl_ann_clust[i].cons_seq[j] + 1;
    if (contig_size == 0)
        return;
    
    // merge close contigs
    vector<bool> temp_reads(recode_data.size(), false);
    vector<int> temp_var(contig_size, 0);
    vector<int64_t> idx;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        if (rl_ann_clust[i].cons_seq.size() == 0)
            continue;
        
        // fill temp_var
        for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k)
            temp_var[this->rl_ann_clust[i].cons_seq[k]] = 1;
        
        // get close neighbor contigs
        vector<int64_t> neighbor_ids;
        unordered_set<int> fp_var;
        unordered_set<int> fn_var;
        unordered_set<int64_t> nn_reads_ids;
        for (auto j = 0; j < this->rl_ann_clust.size(); ++j){
            if (j == i) continue;
            
            // check overlap of two contigs
            int64_t overlap_start = this->rl_ann_clust[i].start >= this->rl_ann_clust[j].start ?
                                    this->rl_ann_clust[i].start : this->rl_ann_clust[j].start;
            int64_t overlap_end = this->rl_ann_clust[i].end <= this->rl_ann_clust[j].end ?
                                  this->rl_ann_clust[i].end : this->rl_ann_clust[j].end;
            int64_t overlap_len = overlap_end - overlap_start + 1;
            if(overlap_len <= 0) continue;
            
            // check similarity
            int64_t n_common = 0;
            int64_t n_union = this->rl_ann_clust[i].cons_seq.size();
            
            for (auto k = 0; k < this->rl_ann_clust[j].cons_seq.size(); ++k){
                if (temp_var[this->rl_ann_clust[j].cons_seq[k]] == 1){
                    ++n_common;
                    temp_var[this->rl_ann_clust[j].cons_seq[k]] = 2;
                }else{
                    ++n_union;
                }
            }
            
            int64_t n_diff = n_union - n_common;
            if (n_diff > max_loci ||
                n_common < int(0.5*(this->rl_ann_clust[i].cons_seq.size()) + 0.5) ||
                n_common < int(0.5*(this->rl_ann_clust[j].cons_seq.size()) + 0.5) ){
                continue;
            }
            
            // record different variants
            neighbor_ids.push_back(j);
            
            for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k){
                if (temp_var[this->rl_ann_clust[i].cons_seq[k]] == 1)
                    fp_var.insert(this->rl_ann_clust[i].cons_seq[k]);
            }
            
            for (auto k = 0; k < this->rl_ann_clust[j].cons_seq.size(); ++k){
                if (temp_var[this->rl_ann_clust[j].cons_seq[k]] == 0)
                    fn_var.insert(this->rl_ann_clust[j].cons_seq[k]);
            }
            
            for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k)
                temp_var[this->rl_ann_clust[i].cons_seq[k]] = 1;
            
            // record neighbor reads
            nn_reads_ids.insert(this->rl_ann_clust[j].nn_reads_id.begin(), this->rl_ann_clust[j].nn_reads_id.end());
        }
        
        // clean temp_var
        for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k)
            temp_var[this->rl_ann_clust[i].cons_seq[k]] = 0;
        
        // quit if no neighbor contigs
        if (neighbor_ids.size() == 0){
            idx.push_back(i);
            for (auto it = nn_reads_ids.begin(); it != nn_reads_ids.end(); ++it)
                temp_reads[*it] = false;
            continue;
        }
        
        // check if the ith contig is from noise
        nn_reads_ids.insert(this->rl_ann_clust[i].nn_reads_id.begin(), this->rl_ann_clust[i].nn_reads_id.end());
        
        for (auto it = nn_reads_ids.begin(); it != nn_reads_ids.end(); ++it)
            temp_reads[*it] = true;

        // test putative false positives
        bool is_noise = false;
        for (auto it = fp_var.begin(); it != fp_var.end(); ++it){
            double cur_count = 0, cur_cvg = 0;
            int64_t cur_locus = int64_t( (*it) / 4);
            if (4*cur_locus+3 >= pu_recode.size() || 4*cur_locus+3 >= pu_recode_ref.size())
                continue;
            
            // get count
            for (auto k = 0; k < pu_recode[*it].size(); ++k){
                if (temp_reads[ pu_recode[*it][k] ])
                    ++cur_count;
            }
            
            // get coverage
            for (int t = 0; t <= 3; ++t){
                for (auto k = 0; k < pu_recode[4*cur_locus+t].size(); ++k){
                    if (temp_reads[ pu_recode[4*cur_locus+t][k] ])
                        ++cur_cvg;
                }
                
                for (auto k = 0; k < pu_recode_ref[4*cur_locus+t].size(); ++k){
                    if (temp_reads[ pu_recode_ref[4*cur_locus+t][k] ])
                        ++cur_cvg;
                }
            }
            
            // test
            if (this->alpha == -1 || this->beta == -1)
                throw runtime_error("this->alpha == -1 || this->beta == -1");
            
            double cur_log_bf = binom_log_bf(cur_count, cur_cvg, this->alpha, this->beta);
            if (cur_log_bf < min_log_bf)
                is_noise = true;
        }
        
        // test putative false negatives
        for (auto it = fn_var.begin(); it != fn_var.end(); ++it){
            double cur_count = 0, cur_cvg = 0;
            int64_t cur_locus = int64_t( (*it) / 4);
            if (4*cur_locus+3 >= pu_recode.size() || 4*cur_locus+3 >= pu_recode_ref.size())
                continue;
            
            // get count and coverage
            for (int t = 0; t <= 3; ++t){
                for (auto k = 0; k < pu_recode[4*cur_locus+t].size(); ++k){
                    if (temp_reads[ pu_recode[4*cur_locus+t][k] ])
                        ++cur_cvg;
                }
                
                for (auto k = 0; k < pu_recode_ref[4*cur_locus+t].size(); ++k){
                    if (temp_reads[ pu_recode_ref[4*cur_locus+t][k] ]){
                        ++cur_count;
                        ++cur_cvg;
                    }
                }
            }
            
            // test
            if (this->alpha == -1 || this->beta == -1)
                throw runtime_error("this->alpha == -1 || this->beta == -1");
            
            double cur_log_bf = binom_log_bf(cur_count, cur_cvg, this->alpha, this->beta);
            if (cur_log_bf < min_log_bf)
                is_noise = true;
        }
        
        if (!is_noise) idx.push_back(i);
        for (auto it = nn_reads_ids.begin(); it != nn_reads_ids.end(); ++it)
            temp_reads[*it] = false;
    }
    
    this->print_rl_ann_clust(out_file, true, idx);
}

void Assembler::correct_contigs(string ann_file, string out_file, double min_overlap_ratio)
{
    this->read_ann_results(ann_file);
    
    // construct hash table fro cons_seq and tested_loci
    vector<unordered_set<int> > cons_seq_hash(this->rl_ann_clust.size(), unordered_set<int>());
    vector<unordered_set<int> > tested_loci_hash(this->rl_ann_clust.size(), unordered_set<int>());
    
    vector<set<int> > cons_seq_corrected(this->rl_ann_clust.size(), set<int>());
    vector<set<int> > tested_loci_corrected(this->rl_ann_clust.size(), set<int>());
    
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        cons_seq_hash[i].insert(this->rl_ann_clust[i].cons_seq.begin(), this->rl_ann_clust[i].cons_seq.end());
        tested_loci_hash[i].insert(this->rl_ann_clust[i].tested_loci.begin(), this->rl_ann_clust[i].tested_loci.end());
    }
    
    // correct contigs
    vector<bool> is_corrected(this->rl_ann_clust.size(), false);
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        
        cons_seq_corrected[i].insert(this->rl_ann_clust[i].cons_seq.begin(), this->rl_ann_clust[i].cons_seq.end());
        tested_loci_corrected[i].insert(this->rl_ann_clust[i].tested_loci.begin(), this->rl_ann_clust[i].tested_loci.end());
        
        for (auto j =0; j < this->rl_ann_clust.size(); ++j){
            if (i == j) continue;
            
            // check overlap ratio
            int overlap_start = this->rl_ann_clust[i].start > this->rl_ann_clust[j].start ? this->rl_ann_clust[i].start : this->rl_ann_clust[j].start;
            int overlap_end = this->rl_ann_clust[i].end < this->rl_ann_clust[j].end ? this->rl_ann_clust[i].end : this->rl_ann_clust[j].end;
            if (overlap_end - overlap_start + 1 < (this->rl_ann_clust[i].end - this->rl_ann_clust[i].start + 1)*min_overlap_ratio )
                continue;
            
            bool is_equal = true;
            // test if cons_seq_j match cons_seq_i
            for (auto k = 0; k < this->rl_ann_clust[j].cons_seq.size(); ++k){
                if ( cons_seq_hash[i].find(this->rl_ann_clust[j].cons_seq[k]) == cons_seq_hash[i].end() &&
                    tested_loci_hash[i].find(this->rl_ann_clust[j].cons_seq[k]/4) != tested_loci_hash[i].end() ){
                    is_equal = false;
                }
            }
            
            // test if cons_seq_i match cons_seq_j
            for (auto k = 0; k < this->rl_ann_clust[i].cons_seq.size(); ++k){
                if ( cons_seq_hash[j].find(this->rl_ann_clust[i].cons_seq[k]) == cons_seq_hash[j].end() ){
                    is_equal = false;
                }
            }
            
            // correct contigs
            if (is_equal){
                is_corrected[i] = true;
                for (auto k = 0; k < this->rl_ann_clust[j].cons_seq.size(); ++k){
                    if (cons_seq_hash[i].find(rl_ann_clust[j].cons_seq[k]) == cons_seq_hash[i].end() &&
                        tested_loci_hash[i].find(rl_ann_clust[j].cons_seq[k]/4) == tested_loci_hash[i].end() &&
                        int(rl_ann_clust[j].cons_seq[k]/4) >= rl_ann_clust[i].start &&
                        int(rl_ann_clust[j].cons_seq[k]/4) <= rl_ann_clust[i].end){
                        
                        cons_seq_corrected[i].insert(rl_ann_clust[j].cons_seq[k]);
                        tested_loci_corrected[i].insert(rl_ann_clust[j].cons_seq[k]/4);
                        
                    }
                }
            }
        }
    }
    
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        if (is_corrected[i]){
            this->rl_ann_clust[i].cons_seq = vector<int>(cons_seq_corrected[i].begin(), cons_seq_corrected[i].end());
            this->rl_ann_clust[i].tested_loci = vector<int>(tested_loci_corrected[i].begin(), tested_loci_corrected[i].end());
        }
    }
    
    this->print_rl_ann_clust(out_file, true);
}



void Assembler::assemble(Graph &gp, string out_ann_file)
{
    // get number of the graph
    size_t nv = get_num_vertices(gp);
    
    if (this->rl_ann_clust.size() == 0){
        this->print_rl_ann_clust(out_ann_file, true);
        return;
    }
        //throw runtime_error("Assembler::assemble(), this->rl_ann_clust.size() == 0");
    
    if (this->rl_ann_clust.size() != nv)
        throw runtime_error("unmatched graph and ann_clust");
    
    // get unambigious path in the graph
    set<GraphPath> upaths = get_unambigious_paths(gp);
    
    // merge consensus sequences in each path
    vector<ConsensusSeq> ann_upath;
    IndexMap index = get(boost::vertex_index, gp);
    for (auto it = upaths.begin(); it != upaths.end(); ++it){
        ConsensusSeq cur_cons;
        cur_cons.seed.push_back(-1);
        cur_cons.neighbors_id.push_back(-1);
        set<int> cur_cons_seq;
        set<int> cur_tested_loci;
        
        for (auto i = 0; i < (*it).size(); ++i){
            int64_t k = index((*it)[i]);
            
            // merge cons_seq
            for (auto j = 0; j < this->rl_ann_clust[k].cons_seq.size(); ++j)
                cur_cons_seq.insert(this->rl_ann_clust[k].cons_seq[j]);
            
            // merge tested_loci
            for (auto j = 0; j < this->rl_ann_clust[k].tested_loci.size(); ++j)
                cur_tested_loci.insert(this->rl_ann_clust[k].tested_loci[j]);
            
            // merge range
            if (i == 0){
                cur_cons.start = this->rl_ann_clust[k].start;
                cur_cons.end = this->rl_ann_clust[k].end;
            }else{
                if (this->rl_ann_clust[k].start < cur_cons.start)
                    cur_cons.start = this->rl_ann_clust[k].start;
                
                if (this->rl_ann_clust[k].end > cur_cons.end)
                    cur_cons.end = this->rl_ann_clust[k].end;
            }
            
        }
        
        for (auto it = cur_cons_seq.begin(); it != cur_cons_seq.end(); ++it)
            cur_cons.cons_seq.push_back(*it);
        
        for (auto it = cur_tested_loci.begin(); it != cur_tested_loci.end(); ++it)
            cur_cons.tested_loci.push_back(*it);
        
        ann_upath.push_back(cur_cons);
    }
    
    this->rl_ann_clust.clear();
    this->rl_ann_clust = ann_upath;
    
    // print ann consensus sequences
    this->print_rl_ann_clust(out_ann_file, true);
    
    // print assembly paths
    ofstream fs_out_ann_file; open_outfile(fs_out_ann_file, out_ann_file + ".upaths");
    for (auto it = upaths.begin(); it != upaths.end(); ++it){
        for (auto i = 0; i < it->size(); ++i){
            fs_out_ann_file << (*it)[i] << ' ';
        }
        fs_out_ann_file << endl;
    }
        
    fs_out_ann_file.close();
}

void Assembler::assemble_unambiguous(IGDA_Graph &gp, string out_file)
{
    // get active vertices (0 in degree)
    stack<int64_t, vector<int64_t> > v_active;
    for (auto it = gp.adj_mat_in.begin(); it != gp.adj_mat_in.end(); ++it){
        if (it->second.size() == 0){
            v_active.push(it->first);
        }
    }
    
    // get unambigious paths (paths where there is no more than two paths connecting its start and end nodes)
    cout << "get unambigious paths (ms)" << endl;
    //set<vector<int64_t> > upaths = get_unambigious_paths_ms(gp);
    set<vector<int64_t> > upaths = get_unambigious_paths_ns(gp);
    
    // merge consensus sequences in each path
    cout << "consensus" << endl;
    vector<ConsensusSeq> ann_upath;
    for (auto it = upaths.begin(); it != upaths.end(); ++it){
        ConsensusSeq cur_cons;
        cur_cons.seed.push_back(-1);
        cur_cons.neighbors_id.push_back(-1);
        set<int> cur_cons_seq;
        set<int> cur_tested_loci;
        
        //for (auto i = 0; i < (*it).size(); ++i){
        for (auto i = 0; i < it->size(); ++i){
            int64_t k = (*it)[i];
            
            // merge cons_seq
            for (auto j = 0; j < this->rl_ann_clust[k].cons_seq.size(); ++j)
                cur_cons_seq.insert(this->rl_ann_clust[k].cons_seq[j]);
            
            // merge tested_loci
            for (auto j = 0; j < this->rl_ann_clust[k].tested_loci.size(); ++j)
                cur_tested_loci.insert(this->rl_ann_clust[k].tested_loci[j]);
            
            // merge range
            if (i == 0){
                cur_cons.start = this->rl_ann_clust[k].start;
                cur_cons.end = this->rl_ann_clust[k].end;
            }else{
                if (this->rl_ann_clust[k].start < cur_cons.start)
                    cur_cons.start = this->rl_ann_clust[k].start;
                
                if (this->rl_ann_clust[k].end > cur_cons.end)
                    cur_cons.end = this->rl_ann_clust[k].end;
            }
            
        }
        
        for (auto it = cur_cons_seq.begin(); it != cur_cons_seq.end(); ++it)
            cur_cons.cons_seq.push_back(*it);
        
        for (auto it = cur_tested_loci.begin(); it != cur_tested_loci.end(); ++it)
            cur_cons.tested_loci.push_back(*it);
        
        ann_upath.push_back(cur_cons);
    }
    this->rl_ann_clust.clear();
    this->rl_ann_clust = ann_upath;
    
    // print ann consensus sequences
    this->print_rl_ann_clust(out_file, true);
    
    // print assembly paths
    ofstream fs_out_ann_file; open_outfile(fs_out_ann_file, out_file + ".upaths");
    for (auto it = upaths.begin(); it != upaths.end(); ++it){
        for (auto i = 0; i < it->size(); ++i){
            fs_out_ann_file << (*it)[i] << ' ';
        }
        fs_out_ann_file << endl;
    }
    
    fs_out_ann_file.close();
}

// polish contigs
void Assembler::polish(string ann_file, string encode_file, string m5_file, string ref_file, string out_file, string tmp_dir, double min_condprob, int min_reads, int min_homo_block_dist, string var_file)
{
    // convert ann consensus sequence to cmpreads
    this->read_ann_results(ann_file);
    ofstream fs_cmpreads_file;
    open_outfile(fs_cmpreads_file, ann_file + ".cmpreads.txt");
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        fs_cmpreads_file << this->rl_ann_clust[i].cons_seq << "," <<endl;
    }
    fs_cmpreads_file.close();
    string cmd = "igda txt2bin " + ann_file + ".cmpreads.txt" + " " + ann_file + ".cmpreads";
    cout << cmd << endl;
    system(cmd.c_str());
    
    // run dforest to calculate conditional probability
    if (var_file ==""){
        cmd = "igda dforest -i -q 2 -r " + to_string(min_reads) + " -f " + to_string(min_condprob) + " -m " + to_string(min_homo_block_dist) + " " +
            encode_file + " " + m5_file + " " + ann_file + ".cmpreads " + ref_file + " " + ann_file + ".dforest " + tmp_dir;
    }else{
        cmd = "igda dforest -i -q 2 -r " + to_string(min_reads) + " -f " + to_string(min_condprob) + " -m " + to_string(min_homo_block_dist) + " -v " + var_file + " " + 
            encode_file + " " + m5_file + " " + ann_file + ".cmpreads " + ref_file + " " + ann_file + ".dforest " + tmp_dir;
    }
    cout << cmd << endl;
    system(cmd.c_str());
    
    /*-------------- polish ann results -------------*/
    // load dforest_data
    unordered_map<int, vector<vector<int>> > dforest_data;
    ifstream fs_dforestfile;
    open_infile(fs_dforestfile, ann_file + ".dforest.all");
    while (true) {
        string buf;
        getline(fs_dforestfile, buf);
        if (fs_dforestfile.eof())
            break;
        
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size() != 7)
            throw runtime_error("permute_encodefile(): buf_vec.size() != 7 in " + ann_file + ".dforest.all");
        dforest_data[stoi(buf_vec[0])].push_back(split_int(buf_vec[6], ','));
    }
    fs_dforestfile.close();
    
    // polish cons_seq of ann
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        if (rl_ann_clust[i].cons_seq.size() <= 1) continue;
        vector<int> cons_seq_polished;
        unordered_set<int> cons_seq_hash(this->rl_ann_clust[i].cons_seq.begin(), this->rl_ann_clust[i].cons_seq.end());
        for (auto j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j){
            auto it = dforest_data.find(rl_ann_clust[i].cons_seq[j]);
            if (it == dforest_data.end()) continue;
            for (auto s = 0; s < it->second.size(); ++s){
                bool is_valid = true;
                for (auto t = 0; t < it->second[s].size(); ++t){
                    auto cur_it = cons_seq_hash.find(it->second[s][t]);
                    if (cur_it == cons_seq_hash.end())
                        is_valid = false;
                }
                if (is_valid){
                    cons_seq_polished.push_back(rl_ann_clust[i].cons_seq[j]);
                    break;
                }
            }
        }
        rl_ann_clust[i].cons_seq = cons_seq_polished;
    }
    
    this->print_rl_ann_clust(out_file, true);
    
}

void Assembler::cut_overhanged_contigs(string ann_file, string ann_cut_file)
{
    this->read_ann_results(ann_file);
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        for (auto j = 0; j < this->rl_ann_clust.size(); ++j){
            if (i == j) continue;
            
            // get overlap region
            int64_t o_start = rl_ann_clust[i].start >= rl_ann_clust[j].start ? rl_ann_clust[i].start : rl_ann_clust[j].start;
            int64_t o_end = rl_ann_clust[i].end <= rl_ann_clust[j].end ? rl_ann_clust[i].end : rl_ann_clust[j].end;
            if (o_start < o_end) continue;
            
            if (this->rl_ann_clust[i].start > this->rl_ann_clust[j].start && this->rl_ann_clust[i].start < this->rl_ann_clust[j].end &&
                this->rl_ann_clust[j].end > this->rl_ann_clust[i].start &&this->rl_ann_clust[j].end < this->rl_ann_clust[i].end){
                
                // check if overlapped region identical
                unordered_set<int> buf;
                for (auto k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k){
                    if (rl_ann_clust[i].cons_seq[k] >= 4*o_start && rl_ann_clust[i].cons_seq[k] <= 4*o_end+3)
                        buf.insert(rl_ann_clust[i].cons_seq[k]);
                }
                
                bool is_identical = true;
                for (auto k = 0; k < rl_ann_clust[j].cons_seq.size(); ++k){
                    if (rl_ann_clust[j].cons_seq[k] >= 4*o_start && rl_ann_clust[j].cons_seq[k] <= 4*o_end+3 &&
                        buf.find(rl_ann_clust[j].cons_seq[k]) == buf.end())
                        is_identical = false;
                }
                
                // if overlapped regions are identical, cut overlaped region for contig i
                if (is_identical){
                    // cut cons_seq
                    vector<int> cons_seq_cut;
                    for (auto k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k){
                        if (rl_ann_clust[i].cons_seq[k] > 4*o_end)
                            cons_seq_cut.push_back(rl_ann_clust[i].cons_seq[k]);
                    }
                    rl_ann_clust[i].cons_seq = cons_seq_cut;
                    
                    // trim tested_loci
                    vector<int> tested_loci_cut;
                    for (auto k = 0; k < rl_ann_clust[i].tested_loci.size(); ++k){
                        if (rl_ann_clust[i].tested_loci[k] > o_end){
                            tested_loci_cut.push_back(rl_ann_clust[i].tested_loci[k]);
                        }
                    }
                    rl_ann_clust[i].tested_loci = tested_loci_cut;
                    
                    // trim start locus
                    rl_ann_clust[i].start = (int)o_end;
                    
                }
            }
        }
    }
    this->print_rl_ann_clust(ann_cut_file);
}

void Assembler::assign_reads_to_contigs(const vector<vector<int> > &recode_data, const vector<ReadRange> &reads_range, bool is_random)
{
    if (recode_data.size() != reads_range.size())
        throw runtime_error("In Assembler::assign_reads_to_contigs, recode_data.size() != reads_range.size()");
    
    size_t genome_size = get_genome_size(reads_range);
    vector<bool> temp_array(genome_size*4+3, false);
    
    ConsensusSeq ref_cons;
    ref_cons.start = 0;
    ref_cons.end = int(genome_size) - 1;
    for (auto i = 0; i < recode_data.size(); ++i){
        // to be removed
        //i = 99774;
        
        if ((i+1) % 1000 == 0)
            cout << i + 1<< endl;
        // get maximal similarity
        vector<double> contig_jaccard;
        vector<double> contig_cvg;
        double max_jaccard = -1;
        for (auto j = 0; j < this->rl_ann_clust.size(); ++j){
            ReadRange contig_range(this->rl_ann_clust[j].start, this->rl_ann_clust[j].end);
            //int min_len = 0.75*(contig_range.second - contig_range.first + 1);
            //double cur_jaccard = sim_jaccard(recode_data[i], this->rl_ann_clust[j].cons_seq, reads_range[i], contig_range, temp_array, true, 0);
            double cur_jaccard = sim_jaccard(recode_data[i], this->rl_ann_clust[j].cons_seq, reads_range[i], contig_range, temp_array, true, 0, true);
            contig_jaccard.push_back(cur_jaccard);
            
            int overlap_start = reads_range[i].first >= contig_range.first ? reads_range[i].first : contig_range.first;
            int overlap_end = reads_range[i].second <= contig_range.second ? reads_range[i].second : contig_range.second;
            int overlap_len = overlap_end - overlap_start + 1;
            
            overlap_len = overlap_len >= 0 ? overlap_len : 0;
            contig_cvg.push_back(overlap_len);
            
            if (cur_jaccard > max_jaccard)
                max_jaccard = cur_jaccard;
        }
        
        if (abs(max_jaccard + 1) <= EPS){
            ref_cons.contig_count += 1.0;
            ref_cons.contig_cvg += reads_range[i].second - reads_range[i].first + 1;
            ref_cons.nn_reads_id.push_back(i);
            continue;
        }
        
        if (max_jaccard <= 0.5)
            continue;
        
        // get contigs having the maximal similarity
        vector<int64_t> contig_id;
        for (auto j = 0; j < this->rl_ann_clust.size(); ++j){
            if ( abs(contig_jaccard[j] - max_jaccard) <= EPS){
                contig_id.push_back(j);
                
            }
        }
        
        for (auto j = 0; j < contig_id.size(); ++j){
            //this->rl_ann_clust[contig_id[j]].contig_count += 1.0/double(contig_id.size());
            //this->rl_ann_clust[contig_id[j]].contig_cvg += contig_cvg[contig_id[j]] / double(contig_id.size());
            this->rl_ann_clust[contig_id[j]].contig_count += 1.0;
            this->rl_ann_clust[contig_id[j]].contig_cvg += contig_cvg[contig_id[j]];
            
            this->rl_ann_clust[contig_id[j]].nn_reads_id.push_back(i);
        }
    }
    
    this->rl_ann_clust.push_back(ref_cons);
    
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i) {
        double cur_contig_len = this->rl_ann_clust[i].end - this->rl_ann_clust[i].start + 1;
        this->rl_ann_clust[i].contig_cvg = this->rl_ann_clust[i].contig_cvg / cur_contig_len;
    }
    
}

void Assembler::find_ncreads(string encode_file, string align_file, string var_file, int topn, double max_dist)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (reads_range.size() != encode_data.size())
        throw runtime_error("reads_range.size() != encode_data.size()");
    
    // get genome size
    size_t genome_size = get_genome_size(reads_range);
    
    // get cumulative sum of variants
    vector<int> var_cdf; get_var_cdf(var_cdf, var_file, genome_size);
    
    // create a template to compare reads
    vector<bool> temp_array(genome_size*4+3, false);
    
    // get topn nearest neighbors for each reads and find non-contained reads (no topn reads can cover all its range)
    nc_reads_id.clear();
    for (auto i = 0; i < encode_data.size(); ++i){
        if ((i+1)%1000 == 0)
            cout << "processed " << i+1 << " / " << encode_data.size() << endl;
        // calculate hamming distance between reads i and other reads
        priority_queue<pair<int,double>, vector<pair<int,double> >, reads_compare_dist > topn_id;
        for (auto j = 0; j < encode_data.size(); ++j){
            if (i == j) continue;
            if (reads_range[i].first >= reads_range[j].second || reads_range[i].second <= reads_range[j].first)
                continue;
            if (reads_range[i].first < reads_range[j].first)
                continue;
            
            double cur_dist = dist_hamming(encode_data[i], encode_data[j], reads_range[i], reads_range[j], var_cdf, temp_array);
            
            if (cur_dist < 0) continue;
            
            topn_id.push(pair<int,double>(j,cur_dist));
        }
        // check if topn reads cover read i
        bool is_nc = true;
        for (auto j = 0; j < topn; ++j){
            if (topn_id.empty()) break;
            
            int cur_id = topn_id.top().first;
            double cur_dist = topn_id.top().second;
            
            if (cur_dist > max_dist)
                break;
            
            if ((reads_range[cur_id].first <= reads_range[i].first && reads_range[cur_id].second > reads_range[i].second) ||
                (reads_range[cur_id].first < reads_range[i].first && reads_range[cur_id].second >= reads_range[i].second) ){
                is_nc = false;
                break;
            }
            topn_id.pop();
        }
        if (is_nc)
            nc_reads_id.push_back(i);
    }
    cout << "processed " << encode_data.size() << " / " << encode_data.size() << endl;

}

bool Assembler::check_pileup(const vector<int> &pu_var_count, const vector<int> &pu_reads_count, int start, int end, const vector<int> &idx, int min_cvg, double min_prop, double max_prop)
{
    if (floor(double(pu_var_count.size()-1) / 4) > (int)pu_reads_count.size() - 1){
        cout << "pu_var_count.size() " << pu_var_count.size() << endl;
        cout << "pu_read_count.size()" << pu_reads_count.size() << endl;
        throw runtime_error("ann_clust: floor(double(pu_var_count.size()-1) / 4) > pu_read_count.size() - 1");
    }
    
    int start_code = 4*start;
    int end_code = 4*end+3;
    end_code = end_code <= pu_var_count.size()-1 ? end_code : (int)pu_var_count.size()-1;
    
    if (start_code >= pu_var_count.size() || end_code >= pu_var_count.size()){
        cout << "start_code = " << start_code << endl;
        cout << "end_code = " << end_code << endl;
        cout << "pu_var_count.size() = " << pu_var_count.size() << endl;
        throw runtime_error("start_code >= pu_var_count.size() || end_code >= pu_var_count.size()");
    }
    
    if (start_code < 0 || end_code < 0){
        cout << "start_code = " << start_code << endl;
        cout << "end_code = " << end_code << endl;
        throw runtime_error("start_code < 0 || end_code < 0");
    }
    
    bool is_homo = true;
    for (auto i = start_code; i < end_code; ++i){
        // calculate variant frequency
        int i_r = int(i/4);
        double cur_prop;
        if (pu_reads_count[i_r] >= min_cvg)
            cur_prop = (double)pu_var_count[i] / pu_reads_count[i_r];
        else
            cur_prop = -1;
        
        // check if current locus is homogeneous
        if (cur_prop > min_prop && cur_prop < max_prop){
            is_homo = false;
            break;
        }
    }
    
    return is_homo;
}

bool Assembler::check_pileup_recode(const vector<int> &pu_var_count, const vector<int> &pu_var_ref_count, int start, int end, const vector<int> &idx, int min_cvg, double min_prop, double max_prop)
{
    // validate inputs
    int start_code = 4*start;
    int end_code = 4*end+3;
    end_code = end_code <= pu_var_count.size()-1 ? end_code : (int)pu_var_count.size()-1;
    
    if (start_code >= pu_var_count.size() || end_code >= pu_var_count.size())
        throw runtime_error("start_code >= pu_var_count.size() || end_code >= pu_var_count.size()");
    
    if (start_code < 0 || end_code < 0)
        throw runtime_error("start_code < 0 || end_code < 0");

    
    bool is_homo = true;
    
    for (auto i = start; i < end; ++i){
        // get coverage of the currecnt locus (only count recoded base including the reference base)
        int64_t cur_cvg = pu_var_count[4*i] + pu_var_count[4*i+1] + pu_var_count[4*i+2] + pu_var_count[4*i+3];
        cur_cvg += pu_var_ref_count[4*i] + pu_var_ref_count[4*i+1] + pu_var_ref_count[4*i+2] + pu_var_ref_count[4*i+3];
        
        if (cur_cvg < min_cvg)
            continue;
        
        // check A, C, G, T
        double cur_prop_A = double(pu_var_count[4*i]) / cur_cvg;
        double cur_prop_C = double(pu_var_count[4*i+1]) / cur_cvg;
        double cur_prop_G = double(pu_var_count[4*i+2]) / cur_cvg;
        double cur_prop_T = double(pu_var_count[4*i+3]) / cur_cvg;
        
        if ((cur_prop_A > min_prop && cur_prop_A < max_prop)||
            (cur_prop_C > min_prop && cur_prop_C < max_prop)||
            (cur_prop_G > min_prop && cur_prop_G < max_prop)||
            (cur_prop_T > min_prop && cur_prop_T < max_prop) ){
            
            is_homo = false;
            break;
        }
    }
   
    return is_homo;
}


void Assembler::find_nccontigs(vector<int64_t> &idx, double min_prop, double min_jaccard, int max_ndiff)
{
    if (this->rl_ann_clust.size() == 0)
        return;
    
    // init jaccard
    size_t gsize = 0;
    for (auto i = 0; i < this->rl_ann_clust.size(); ++i){
        if (this->rl_ann_clust[i].end > gsize)
            gsize = this->rl_ann_clust[i].end;
    }
    vector<bool> temp_array(gsize*4 + 3, false);
    
    // get maximal encoded cons_seq and tested_loci
    int64_t temp_size = 0;
    int64_t temp_tested_size = 0;
    for (int64_t i = 0; i < this->rl_ann_clust.size(); ++i){
        for (int64_t j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j)
            if (rl_ann_clust[i].cons_seq[j] > temp_size)
                temp_size = rl_ann_clust[i].cons_seq[j];
        
        for (int64_t j = 0; j < this->rl_ann_clust[i].tested_loci.size(); ++j)
            if (rl_ann_clust[i].tested_loci[j] > temp_tested_size)
                temp_tested_size = rl_ann_clust[i].tested_loci[j];
        
    }
    if (temp_size == 0)
        throw runtime_error("Assembler::find_nccontigs, temp_size == 0");
    if (temp_tested_size == 0)
        throw runtime_error("Assembler::find_nccontigs, temp_tested_size == 0");
    
    ++temp_size;
    ++temp_tested_size;
    
    // generate a template vector, note that temp_vec is as least 4 times larger than temp_tested_vec
    vector<bool> temp_vec_i(temp_size, false);
    vector<bool> temp_vec_j(temp_size, false);
    vector<bool> temp_tested_vec_i(temp_tested_size, false);
    vector<bool> temp_tested_vec_j(temp_tested_size, false);
    
    // pairwise compare contigs
    for (int64_t i = 0; i < rl_ann_clust.size(); ++i){
        // fill in template by the ith cons_seq and tested_loci
        for (int64_t k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k)
            temp_vec_i[rl_ann_clust[i].cons_seq[k]] = true;
        
        for (int64_t k = 0; k < rl_ann_clust[i].tested_loci.size(); ++k)
            temp_tested_vec_i[rl_ann_clust[i].tested_loci[k]] = true;
        
        // compare contigs
        bool is_nc = true;
        int start_code = 4*rl_ann_clust[i].start;
        int end_code = 4*rl_ann_clust[i].end +3;
        
        for (int64_t j = 0; j < this->rl_ann_clust.size(); ++j){
            if (j == i) continue;
            
            if (rl_ann_clust[i].start == rl_ann_clust[j].start && rl_ann_clust[i].end == rl_ann_clust[j].end && j < i)
                continue;
            
            if (rl_ann_clust[j].start <= rl_ann_clust[i].start && rl_ann_clust[j].end >= rl_ann_clust[i].end){
            //if ((rl_ann_clust[j].start <= rl_ann_clust[i].start && rl_ann_clust[j].end > rl_ann_clust[i].end) || (rl_ann_clust[j].start < rl_ann_clust[i].start && rl_ann_clust[j].end >= rl_ann_clust[i].end) ){
                // calculate jaccard index
                double cur_jaccard = sim_jaccard(rl_ann_clust[i].cons_seq, rl_ann_clust[j].cons_seq,
                                                 ReadRange(rl_ann_clust[i].start, rl_ann_clust[i].end),
                                                 ReadRange(rl_ann_clust[j].start, rl_ann_clust[j].end),
                                                 temp_array);
                
                double n_cons_seq_i = rl_ann_clust[i].cons_seq.size();
                double n_cons_seq_j = 0;
                double n_overlap = 0;
                // fill in template by the jth cons_seq and tested_loci
                bool is_diff = false;
                int64_t n_diff = 0;
                for (int64_t k = 0; k < rl_ann_clust[j].cons_seq.size(); ++k){
                    temp_vec_j[rl_ann_clust[j].cons_seq[k]] = true;
                    if (!temp_vec_i[rl_ann_clust[j].cons_seq[k]] && temp_tested_vec_i[rl_ann_clust[j].cons_seq[k] / 4]){
                        is_diff = true;
                        ++n_diff;
                    }
                    
                    if (temp_vec_i[rl_ann_clust[j].cons_seq[k]])
                        ++n_overlap;
                    
                    if (rl_ann_clust[j].cons_seq[k] >= start_code && rl_ann_clust[j].cons_seq[k] <= end_code)
                        ++n_cons_seq_j;
                }
                
                for (int64_t k = 0; k < rl_ann_clust[j].tested_loci.size(); ++k){
                    temp_tested_vec_j[rl_ann_clust[j].tested_loci[k]] = true;
                }
                
                // scan the ith contig again
                for (int64_t k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k){
                    if (!temp_vec_j[rl_ann_clust[i].cons_seq[k]] && temp_tested_vec_j[rl_ann_clust[i].cons_seq[k] / 4]){
                        is_diff = true;
                        ++n_diff;
                    }
                }
                
                //if ((!is_diff || (cur_jaccard >= min_jaccard && n_diff <= 1)) && n_overlap > min_prop*n_cons_seq_i && n_overlap > min_prop*n_cons_seq_j)
                if ((!is_diff || (cur_jaccard >= min_jaccard && n_diff <= max_ndiff)) && n_overlap >= min_prop*n_cons_seq_i && n_overlap >= min_prop*n_cons_seq_j)
                //if ((!is_diff || (cur_jaccard >= min_jaccard || n_diff <= 1)) && n_overlap > min_prop*n_cons_seq_i && n_overlap > min_prop*n_cons_seq_j)
                    is_nc = false;
                
                // clear template of the jth contig
                for (int64_t k = 0; k < rl_ann_clust[j].cons_seq.size(); ++k){
                    temp_vec_j[rl_ann_clust[j].cons_seq[k]] = false;
                }
                
                for (int64_t k = 0; k < rl_ann_clust[j].tested_loci.size(); ++k){
                    temp_tested_vec_j[rl_ann_clust[j].tested_loci[k]] = false;
                }
                
            }
            if(!is_nc) break;
        }
        if (is_nc)
            idx.push_back(i);
        
        
        // clear template
        for (int64_t k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k)
            temp_vec_i[rl_ann_clust[i].cons_seq[k]] = false;
        
        for (int64_t k = 0; k < rl_ann_clust[i].tested_loci.size(); ++k)
            temp_tested_vec_i[rl_ann_clust[i].tested_loci[k]] = false;
        
    }
}

void Assembler::find_nccontigs_legacy(vector<int64_t> &idx)
{
    // get maximal encoded cons_seq
    int64_t temp_size = 0;
    for (int64_t i = 0; i < this->rl_ann_clust.size(); ++i)
        for (int64_t j = 0; j < this->rl_ann_clust[i].cons_seq.size(); ++j)
            if (rl_ann_clust[i].cons_seq[j] > temp_size)
                temp_size = rl_ann_clust[i].cons_seq[j];
    ++temp_size;
    
    // generate a template vector
    vector<bool> temp_vec(temp_size, false);
    
    // pairwise compare cons_seq
    for (int64_t i = 0; i < this->rl_ann_clust.size(); ++i){
        // fill in template by the ith cons_seq
        for (auto k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k)
            temp_vec[rl_ann_clust[i].cons_seq[k]] = true;
        
        bool is_nc = true;
        int start_code = 4*rl_ann_clust[i].start;
        int end_code = 4*rl_ann_clust[i].end +3;
        for (int64_t j = 0; j < this->rl_ann_clust.size(); ++j){
            if (j == i) continue;
            if ((rl_ann_clust[j].start <= rl_ann_clust[i].start && rl_ann_clust[j].end > rl_ann_clust[i].end) || (rl_ann_clust[j].start < rl_ann_clust[i].start && rl_ann_clust[j].end >= rl_ann_clust[i].end) ){
                int n_match = 0;
                for (auto k = 0; k < rl_ann_clust[j].cons_seq.size(); ++k){
                    if (rl_ann_clust[j].cons_seq[k] >= start_code && rl_ann_clust[j].cons_seq[k] <= end_code){
                        if (temp_vec[rl_ann_clust[j].cons_seq[k]])
                            ++n_match;
                        else
                            break;
                    }
                }
                if (n_match == rl_ann_clust[i].cons_seq.size())
                    is_nc = false;
            }
            if(!is_nc) break;
        }
        if (is_nc)
            idx.push_back(i);
        
        // clear template
        for (auto k = 0; k < rl_ann_clust[i].cons_seq.size(); ++k)
            temp_vec[rl_ann_clust[i].cons_seq[k]] = false;
    }
}

void Assembler::get_consensus_recode(ConsensusSeq &cons, const vector<int> &pu_var_count, const vector<int> &pu_var_ref_count, int start, int end, int min_cvg)
{
    cons.cons_seq.clear();
    cons.tested_loci.clear();
    cons.start = start;
    cons.end = start;
    vector<double> prop(pu_var_count.size(),-1);
    
    // validate inputs
    if (pu_var_count.size()==0)
        return;
    
    int start_code = 4*start;
    int end_code = 4*end+3;
    end_code = end_code <= pu_var_count.size()-1 ? end_code : (int)pu_var_count.size()-1;
    
    if (start_code >= pu_var_count.size() || end_code >= pu_var_count.size())
        throw runtime_error("start_code >= pu_var_count.size() || end_code >= pu_var_count.size()");
    
    if (start_code < 0 || end_code < 0)
        throw runtime_error("start_code < 0 || end_code < 0");
    
    for (auto i = start; i < end; ++i){
        // get coverage of the currecnt locus (only count recoded base including the reference base)
        int64_t cur_cvg = pu_var_count[4*i] + pu_var_count[4*i+1] + pu_var_count[4*i+2] + pu_var_count[4*i+3];
        cur_cvg += pu_var_ref_count[4*i] + pu_var_ref_count[4*i+1] + pu_var_ref_count[4*i+2] + pu_var_ref_count[4*i+3];
        
        if (cur_cvg < min_cvg)
            continue;
        
        cons.end = i;
        cons.tested_loci.push_back(i);
        
        // check A, C, G, T
        double cur_prop_A = double(pu_var_count[4*i]) / cur_cvg;
        double cur_prop_C = double(pu_var_count[4*i+1]) / cur_cvg;
        double cur_prop_G = double(pu_var_count[4*i+2]) / cur_cvg;
        double cur_prop_T = double(pu_var_count[4*i+3]) / cur_cvg;
        
        if (cur_prop_A > 0.5)
            cons.cons_seq.push_back(4*i); 
        
        if (cur_prop_C > 0.5)
            cons.cons_seq.push_back(4*i+1);
        
        if (cur_prop_G > 0.5)
            cons.cons_seq.push_back(4*i+2);
        
        if (cur_prop_T > 0.5)
            cons.cons_seq.push_back(4*i+3);
    }
    
}

/*void Assembler::correct_contigs(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range, const vector<int> &var_cdf, const vector<bool> &temp_array, int min_cvg, double min_prop, double max_prop)
{
    vector<bool> temp_array_dual(temp_array.size(), false);
    
    // correct each contig
    for (auto i = 0; i < rl_ann_clust.size(); ++i){
        if (rl_ann_clust[i].cons_seq.size() == 0)
            continue;
        
        for (auto j = 0; j < rl_ann_clust.size(); ++j){
            if (j == i) continue;
        }
    }
}*/


void Assembler::print_correct_reads_raw(const CmpreadsDiffRead &cmpread, ofstream &fs_testfile)
{
    for (int i = 0; i < (int)cmpread.cmpreads_diff.size(); ++i){
        fs_testfile << cmpread.read_id << '\t' << cmpread.cmpreads_diff[i].start << '\t' << cmpread.cmpreads_diff[i].end << '\t';
        fs_testfile << cmpread.cmpreads_diff[i].cand_loci << '\t' << cmpread.cmpreads_diff[i].cand_loci_diff << '\t';
        fs_testfile << cmpread.cmpreads_diff[i].condprob << '\t' << cmpread.cmpreads_diff[i].condprob_diff << endl;
    }
}

void Assembler::print_correct_reads(const CmpreadsDiffRead &cmpread, ofstream &fs_outfile)
{
    for (auto it = cmpread.encode_corrected.begin(); it != cmpread.encode_corrected.end(); ++it){
        fs_outfile << *it << '\t';
    }
    fs_outfile << endl;
}




