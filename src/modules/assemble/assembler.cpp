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
            
            if (jaccard_index < min_jaccard_index)
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


void Assembler::ref_reconstruct(const stxxl::vector<Align> &align_data, string &ref_name, string &ref_seq)
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

void Assembler::ann_clust_recode(string recode_file, string recode_ref_file, string align_file, string var_file, int min_cvg, double min_prop, double max_prop, int topn, int max_nn, double max_dist)
{
    /*------------ find nc-reads -----------*/
    cout << "find non-contained reads" << endl;
    this->find_ncreads(recode_file, align_file, var_file, topn, max_dist);
    //this->nc_reads_id = {17038};
    cout << "number of nc-reads: " << nc_reads_id.size() << endl;
    
    /*------------ use nc-reads seed to cluster ----------*/
    cout << "use non-contained reads as seed to cluster" << endl;
    
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
    
    if (reads_range.size() != recode_data.size())
        throw runtime_error("reads_range.size() != recode_data.size()");
    
    if (reads_range.size() != recode_ref_data.size())
        throw runtime_error("reads_range.size() != recode_ref_data.size()");
    
    // get genome size
    size_t genome_size = get_genome_size(reads_range);
    
    // get cumulative sum of variants
    vector<int> var_cdf; get_var_cdf(var_cdf, var_file, genome_size);
    
    // create a template to compare reads
    vector<bool> temp_array(genome_size*4+3, false);
    
    vector<int> cur_pu_var_count(4*(genome_size-1)+3+1, 0);
    vector<int> cur_pu_var_ref_count(4*(genome_size-1)+3+1, 0);
    vector<int> cur_pu_reads_count(genome_size, 0);
    
    // get topn nearest neighbors for each reads and find non-contained reads (no topn reads can cover all its range)
    int64_t n_nc_reads = 0;
    for (auto i : nc_reads_id){
        ++n_nc_reads;
        if (n_nc_reads % 1000 == 0)
            cout << "processed " << n_nc_reads << " / " << nc_reads_id.size() << endl;
        
        // calculate hamming distance between reads i and other reads and get topn nearest neighbors
        //priority_queue<pair<int,double>, vector<pair<int,double> >, reads_compare_dist > topn_id;
        priority_queue<pair<int,double>, vector<pair<int,double> >, reads_compare_sim > topn_id;
        for (auto j = 0; j < recode_data.size(); ++j){
            if (i == j) continue;
            if (reads_range[i].first >= reads_range[j].second || reads_range[i].second <= reads_range[j].first)
                continue;
            if (reads_range[i].first < reads_range[j].first)
                continue;
            
            //double cur_dist = dist_hamming(recode_data[i], recode_data[j], reads_range[i], reads_range[j], var_cdf, temp_array);
            double cur_dist = sim_jaccard(recode_data[i], recode_data[j], reads_range[i], reads_range[j], temp_array, true);
            
            if (cur_dist < 0) continue;
            
            topn_id.push(pair<int,double>(j,cur_dist));
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
        
        // to be removed
        cout << cur_neighbors_topn << endl;
        cout << cur_neighbors << endl;
        
        // pileup all neighbors and pop from most distant neighbor until all loci are homogeneous
        // pileup topn neighbors
        
        unordered_set<int64_t> mod_idx_var;
        unordered_set<int64_t> mod_idx_var_ref;
        for (auto j : cur_neighbors){
            pileup_var_online_count(cur_pu_var_count, recode_data[j], mod_idx_var);
            pileup_var_online_count(cur_pu_var_ref_count, recode_ref_data[j], mod_idx_var_ref);
        }
        
        // check pileup of all the neighbors
        bool is_homo = this->check_pileup_recode(cur_pu_var_count, cur_pu_var_ref_count, reads_range[i].first, reads_range[i].second, vector<int>(), min_cvg, min_prop, max_prop);
        
        // if not all loci are homogeneous, pop neighbors from the most distant one until all loci are homogeneous or number of neighbors <= topn
        if (!is_homo){
            for (auto j = 0; j < cur_neighbors.size(); ++j){
                int t = (int)cur_neighbors.size() - 1 - j;
                if (t < topn) break;
                pileup_var_online_count_pop(cur_pu_var_count, recode_data[cur_neighbors[t]]);
                pileup_var_online_count_pop(cur_pu_var_ref_count, recode_ref_data[cur_neighbors[t]]);
                is_homo = this->check_pileup_recode(cur_pu_var_count, cur_pu_var_ref_count, reads_range[i].first, reads_range[i].second, vector<int>(), min_cvg, min_prop, max_prop);
                if (is_homo) break;
            }
        }
        
        if (is_homo){
            ConsensusSeq cur_cons;
            cur_cons.seed = recode_data[i];
            cur_cons.neighbors_id = cur_neighbors;
            this->get_consensus_recode(cur_cons, cur_pu_var_count, cur_pu_var_ref_count, reads_range[i].first, reads_range[i].second, min_cvg);
            rl_ann_clust.push_back(cur_cons);
        }
        
        
        // clean cur_pu_var_count and cur_pu_reads_count
        for (auto it = mod_idx_var.begin(); it != mod_idx_var.end(); ++it)
            cur_pu_var_count[*it] = 0;
        for (auto it = mod_idx_var_ref.begin(); it != mod_idx_var_ref.end(); ++it)
            cur_pu_var_ref_count[*it] = 0;
        
    }
    
    
}

void ann_clust_recode_recursive(string recode_file, string recode_ref_file, string align_file, string var_file, int min_cvg = 20, double min_prop = 0.2, double max_prop = 0.7, int topn = 30, int max_nn = 200, int max_iter = 5)
{
    
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
        fs_outfile << this->rl_ann_clust[i].cons_seq << '\t';
        fs_outfile << this->rl_ann_clust[i].start << '\t' << this->rl_ann_clust[i].end;
        
        if (is_metric){
            fs_outfile << '\t';
            
            fs_outfile << this->rl_ann_clust[i].seed << '\t';
            
            /*int start_code = 4*rl_ann_clust[i].start;
            int end_code = 4*rl_ann_clust[i].end + 3;
            end_code = end_code <= (int)rl_ann_clust[i].pu_var_count.size()-1?  end_code : (int)rl_ann_clust[i].pu_var_count.size()-1;
            
            for (auto j = start_code; j < end_code; ++j)
                if (rl_ann_clust[i].pu_var_count[j] > 0)
                    fs_outfile << j << ',';
            fs_outfile << '\t';
            
            for (auto j = start_code; j < end_code; ++j)
                if (rl_ann_clust[i].pu_var_count[j] > 0)
                    fs_outfile << rl_ann_clust[i].prop[j] << ',';
            fs_outfile << '\t';
            
            for (auto j = start_code; j < end_code; ++j)
                if (rl_ann_clust[i].pu_var_count[j] > 0)
                    fs_outfile << rl_ann_clust[i].pu_var_count[j] << ',';
            fs_outfile << '\t';
            
            for (auto j = start_code; j < end_code; ++j){
                if (rl_ann_clust[i].pu_var_count[j] > 0){
                    int j_r = int(j/4);
                    fs_outfile << rl_ann_clust[i].pu_read_count[j_r] << ',';
                }
            }
            fs_outfile << '\t';*/
            
            fs_outfile << rl_ann_clust[i].neighbors_id;
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


void Assembler::find_nccontigs(vector<int64_t> &idx)
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




