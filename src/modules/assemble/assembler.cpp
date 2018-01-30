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
    loadreadsrange(reads_range, align_file, 'm');
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
    loadreadsrange(reads_range, align_file, 'm');
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
        ++cur_count;
        for (int i=pre_pos; i<=cur_pos; ++i)
            cdf_var[i] = cur_count;
        pre_pos = cur_pos + 1;
    }
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



void Assembler::jaccard_index(string encode_file, string align_file, string out_file)
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
            
            fs_out << i <<',' << j << ',' << jaccard_index << ',';
            fs_out << n_intersect << ',' << n_union << ',' << n_overlap << ',';
            fs_out << start << ',' << end << ',' << code_start << ',' << code_end << endl;
            
            fs_out << j <<',' << i << ',' << jaccard_index << ',';
            fs_out << n_intersect << ',' << n_union << ',' << n_overlap << ',';
            fs_out << start << ',' << end << ',' << code_start << ',' << code_end << endl;
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();
}

void Assembler::assemble(string encode_file, string align_file, vector<vector<int> > &centroid, vector<ReadRange> &centroid_range,
                         vector<vector<int> > &idx_on, vector<int> &n_idx_on, int min_idx_on, int min_overlap, int max_iter)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load read range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    // load alignment data
    AlignReaderM5 AlignReaderM5_obj;
    stxxl::vector<Align> align_data;
    AlignReaderM5_obj.read(align_file, align_data);
    
    // reconstruct reference genome from align_data
    string ref_name; string ref_seq;
    this->ref_reconstruct(align_data, ref_name, ref_seq);
    
    /*----------- first round of matrix decomposition -----------*/
    this->assemble_core(encode_data, reads_range, centroid, centroid_range, idx_on, n_idx_on, min_idx_on, min_overlap, max_iter);
    
    /*----------- realign each centroid--------*/
    Alignment aligner;
    AlignCoderSNV aligncoder;
    for (int i=82; i<(int)centroid.size(); ++i){
        cout << i << endl;
        // construct haplotype sequence
        string haplo_seq;
        this->haplo_seq_construct(centroid[i], ref_seq, haplo_seq);
        
        // realign each read (belongs to the current haplotype) to the haplotype sequence
        vector<vector<int> > haplo_encode_data;
        vector<ReadRange> haplo_reads_range;
        StripedSmithWaterman::Alignment result;
        for (int j = 0; j<(int)idx_on[i].size(); j++){
            string cur_qSeq = align_data[idx_on[i][j]].qSeq;
            string cur_tSeq = haplo_seq.substr(align_data[idx_on[i][j]].tStart, align_data[idx_on[i][j]].tEnd - align_data[idx_on[i][j]].tStart + 1);
            
            // align
            aligner.local_align(cur_qSeq, cur_tSeq, result);
           
            // encode
            vector<int> cur_encode_data;
            aligncoder.encode(result, cur_qSeq, cur_tSeq, align_data[i].tStart + result.ref_begin, cur_encode_data);
            haplo_encode_data.push_back(cur_encode_data);
            
            // add realignment
        }
    }
    
}


void Assembler::assemble_core(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                              vector<vector<int> > &centroid, vector<ReadRange> &centroid_range, vector<vector<int> > &idx_on,
                              vector<int> &n_idx_on, int min_idx_on, int min_overlap, int max_iter)
{
    // get non-contained reads
    vector<bool> is_contained = this->check_contained_reads(encode_data, reads_range, min_overlap, true);
    vector<vector<int> > centroid_seed;
    vector<ReadRange> centroid_range_seed;
    for (int i=0; i<(int)is_contained.size(); ++i){
        if (!is_contained[i]){
            centroid_seed.push_back(encode_data[i]);
            centroid_range_seed.push_back(reads_range[i]);
        }
    }
    
    // rank 1 matrix facterization for each seed
    
    for (int i=0; i<(int)centroid_seed.size(); ++i){
        vector<int> cur_centroid = centroid_seed[i];
        ReadRange cur_centroid_range = centroid_range_seed[i];
        vector<int> cur_idx_on; vector<int> idx_off;
        this->mat_fac_rank_1(encode_data, reads_range, cur_centroid, cur_centroid_range, cur_idx_on, idx_off, min_idx_on, min_overlap, max_iter);
        centroid.push_back(cur_centroid);
        centroid_range.push_back(cur_centroid_range);
        n_idx_on.push_back((int)cur_idx_on.size());
        idx_on.push_back(cur_idx_on);
    }
    
}


vector<bool> Assembler::check_contained_reads(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                                   int min_overlap, bool rm_empty_centroid)
{
    if (encode_data.size() != reads_range.size())
        throw runtime_error("encode_data.size() != reads_range.size()");
    
    // setup template and counter
    int temp_size = 0;
    for (int i=0; i<(int)encode_data.size(); ++i)
        for (int j=0; j<(int)encode_data[i].size(); ++j)
            if (encode_data[i][j] > temp_size)
                temp_size = encode_data[i][j];
    temp_size = temp_size + 1;
    vector<int> temp_vec(temp_size, 0);
    int counter = 0;
    
    // check each read to see if they are contained
    vector<bool> is_contained(encode_data.size(), false);
    for (int i=0; i<(int)encode_data.size(); ++i){
        if (encode_data[i].size()==0){
            if (rm_empty_centroid)
                is_contained[i] = true;
            continue;
        }

        // make sure counter does not exceed
        if (counter >= numeric_limits<int>::max()-1)
            throw runtime_error("counter >= numeric_limits<int>::max()-1");
        
        // fill in temp_vec with encode_data[i]
        for (int j=0; j<(int)encode_data[i].size(); ++j)
            temp_vec[encode_data[i][j]] = counter;
        
        // check overlap between encode_data[i] and encode_data[j]
        for (int j=0; j<(int)encode_data.size(); ++j){
            if (j == i)
                continue;
            
            // check if range of read i is contained in range of read j
            if (reads_range[i].first < reads_range[j].first || reads_range[i].second > reads_range[j].second)
                continue;
            
            // calculate overlap
            int overlap_start = reads_range[i].first;
            int overlap_end = reads_range[i].second;
            int overlap_len = overlap_end - overlap_start + 1;
            if (overlap_len < min_overlap)
                continue;
            
            // calculate number of variants of encode_data[j] overlapping with the reads[i]
            int cur_reads_range_code_start = 4*overlap_start;
            int cur_reads_range_code_end = 4*overlap_end + 3;
            int n_overlap = 0;
            for (int k=0; k<(int)encode_data[j].size(); ++k)
                if (encode_data[j][k] >= cur_reads_range_code_start && encode_data[j][k] <= cur_reads_range_code_end)
                    ++n_overlap;
            if (n_overlap == 0)
                continue;
            
            // calculate number of matches
            int n_match = 0;
            for (int k=0; k<(int)encode_data[j].size(); ++k)
                if (temp_vec[encode_data[j][k]] == counter)
                    ++n_match;
            
            // check if read i is contained in read j
            if (n_match>= ceil(double(n_overlap) / 2)){
                is_contained[i] = true;
                break;
            }
            

        }
        
        ++counter;
    }
    
    return is_contained;
}

int Assembler::mat_fac_rank_1(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                    vector<int> &centroid, const ReadRange &centroid_range,
                    vector<int> &idx_on, vector<int> &idx_off, int min_idx_on, int min_overlap, int max_iter)
{
    if (encode_data.size() != reads_range.size())
        throw runtime_error("encode_data.size() != reads_range.size()");
    
    if (centroid.size()==0)
        return 0;
    
    vector<int> old_centroid = centroid;
    vector<int> old_idx_on = idx_on;
    vector<int> old_idx_off = idx_off;
    
    int n_iter = 0;
    for (int i=0; i<=max_iter; ++i){
        vector<int> new_centroid = old_centroid;
        vector<int> new_idx_on, new_idx_off;
        mat_fac_rank_1_core(encode_data, reads_range, new_centroid, centroid_range,
                            new_idx_on, new_idx_off, min_overlap);
        
        ++n_iter;

        // stop if centroid has no change or too few reads match
        if (new_centroid == old_centroid || new_idx_on.size() < min_idx_on){
            old_idx_on = new_idx_on;
            old_idx_off = new_idx_off;
            break;
        }
        
        // if not converge yet, update old_centroid and go on
        old_centroid = new_centroid;
        old_idx_on = new_idx_on;
        old_idx_off = new_idx_off;
    }
    
    centroid = old_centroid;
    idx_on = old_idx_on;
    idx_off = old_idx_off;
    
    return n_iter;
}

void Assembler::mat_fac_rank_1_core(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                                    vector<int> &centroid, const ReadRange &centroid_range,
                                    vector<int> &idx_on, vector<int> &idx_off,  int min_overlap)
{
    // setup template vector for centroid
    int temp_size = *max_element(begin(centroid), end(centroid)) + 1;
    vector<bool> temp_vec(temp_size, false);
    for (int i=0; i<(int)centroid.size(); ++i)
        temp_vec[centroid[i]] = true;
    
    // match all the reads to the template vector
    for (int i=0; i<(int)encode_data.size(); ++i){
        // calculate overlap
        int overlap_start = centroid_range.first >= reads_range[i].first ? centroid_range.first : reads_range[i].first;
        int overlap_end = centroid_range.second <= reads_range[i].second ? centroid_range.second : reads_range[i].second;
        int overlap_len = overlap_end - overlap_start + 1;
        if (overlap_len < min_overlap)
            continue;
        
        // calculate number of variants of centroid overlapping with the current reads
        int cur_reads_range_code_start = 4*reads_range[i].first;
        int cur_reads_range_code_end = 4*reads_range[i].second + 3;
        int n_overlap = 0;
        for (int j=0; j<(int)centroid.size(); ++j)
            if(centroid[j]>=cur_reads_range_code_start && centroid[j]<=cur_reads_range_code_end)
                n_overlap++;
        if (n_overlap == 0)
            continue;
        
        // calculate number of matches
        int n_match = 0;
        for (int j=0; j<(int)encode_data[i].size(); ++j){
            if (encode_data[i][j]>=temp_size)
                break;
            if (temp_vec[encode_data[i][j]])
                ++n_match;
        }
        if (n_match>centroid.size())
            throw runtime_error("n_match > centroid.size() at line " + to_string(i));
        
        // if #shared variants between reads and centroid >= 50%, then they are grouped and index is put in idx_on otherwise idx_off
        if (n_match>= ceil(double(n_overlap) / 2))
            idx_on.push_back(i);
        
    }
    
    // recalculate new centroid
    int centroid_range_code_start = 4*centroid_range.first;
    int centroid_range_code_end = 4*centroid_range.second + 3;
    temp_size = 4*centroid_range.second+3 + 1;
    vector<int> temp_vec_var(temp_size, 0);
    vector<int> temp_vec_reads(temp_size, 0);
    set<int> var_list;
    
    // pileup variants and count number of variants in each position
    for (int i=0; i<(int)idx_on.size(); ++i){
        
        for (int j=0; j<(int)encode_data[idx_on[i]].size(); ++j){
            if (encode_data[idx_on[i]][j] < centroid_range_code_start)
                continue;
            if (encode_data[idx_on[i]][j] > centroid_range_code_end)
                break;
            ++temp_vec_var[ encode_data[idx_on[i]][j] ];
            var_list.insert(encode_data[idx_on[i]][j]);
        }
    }
    
    // for each variants in var_list, cacculate number of reads covering it.
    vector<int> var_list_vec(var_list.size(), 0);
    int vec_i = 0;
    for (auto it=var_list.begin(); it!=var_list.end(); ++it){
        var_list_vec[vec_i] = *it; ++vec_i;
    }
    
    for (int i=0; i<(int)idx_on.size(); ++i){
        int cur_reads_range_code_start = 4*reads_range[idx_on[i]].first;
        int cur_reads_range_code_end = 4*reads_range[idx_on[i]].second + 3;
        for (int j=0; j<(int)var_list_vec.size(); ++j){
            if (var_list_vec[j] >= cur_reads_range_code_start &&
                var_list_vec[j] <= cur_reads_range_code_end)
                ++temp_vec_reads[ var_list_vec[j] ];
        }
    }
    
    // calculate new centroid
    vector<int> new_centroid;
    for (int i=0; i<(int)var_list_vec.size(); ++i){
        if (temp_vec_var[ var_list_vec[i] ] >= ceil(double(temp_vec_reads[var_list_vec[i]]) / 2))
            new_centroid.push_back(var_list_vec[i]);
    }
   
    centroid = new_centroid;
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




