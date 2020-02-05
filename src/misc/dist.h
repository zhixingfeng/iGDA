//
//  dist.h
//  iGDA
//
//  Created by Zhixing Feng on 6/15/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#ifndef iGDA_dist_h
#define iGDA_dist_h

inline double sim_jaccard(const vector<int> &encode_1, const vector<int> &encode_2, const ReadRange &range_1, const ReadRange &range_2, vector<bool> &temp_array, bool check_ref = true, int min_overlap = 500, bool is_asym = false)
{
    ReadRange range_overlap(range_1.first >= range_2.first ? range_1.first : range_2.first, range_1.second <= range_2.second ? range_1.second : range_2.second);
    
    // get overlap length
    int overlap = range_overlap.second - range_overlap.first + 1;
    if (overlap < min_overlap)
        return -1;
    
    int n_intersect = 0;
    int n_union = 0;
    int n_overlap_1 = 0;
    int n_overlap_2 = 0;
    // read 1
    for (auto i = 0; i < encode_1.size(); ++i){
        if (is_asym){
            temp_array[encode_1[i]] = true;
            ++n_union;
            ++n_overlap_1;
        }else{
            if (encode_1[i] >= 4*range_overlap.first && encode_1[i] <= 4*range_overlap.second+3){
                temp_array[encode_1[i]] = true;
                ++n_union;
                ++n_overlap_1;
            }
        }
    }
    
    if (n_union == 0)
        return -1;
    
    // read_2
    for (auto i = 0; i < encode_2.size(); ++i){
        if (encode_2[i] >= 4*range_overlap.first && encode_2[i] <= 4*range_overlap.second+3){
            if (temp_array[encode_2[i]])
                ++n_intersect;
            else
                ++n_union;
            ++n_overlap_2;
        }
    }
    
    // clear up temp_array
    for (auto i = 0; i < encode_1.size(); ++i)
        temp_array[encode_1[i]] = false;
    
    if (check_ref){
        if (2*n_intersect <= encode_1.size())
            return -1;
        /*if (2*n_intersect < n_overlap_1 || 2*n_intersect < n_overlap_2){
            return -1;
        }*/
    }
    return (double)n_intersect / n_union;
}

// mutual information weighted hamming distance
inline void cal_locus_specific_mi(const vector<int> &pu_var_count, const vector<int> &pu_var_ref_count, const vector<int> &var_encode,
                                  vector<double> &weights_11, vector<double> &weights_10, vector<double> &weights_01, vector<double> &weights_00)
{
    weights_11.clear();
    weights_10.clear();
    weights_01.clear();
    weights_00.clear();
    
    
}

inline double dist_hamming_mi(const vector<int> &encode_1, const vector<int> &encode_2, const ReadRange &range_1, const ReadRange &range_2,
                           const vector<int> &var_cdf, vector<bool> &temp_array, int min_overlap = 500)
{
    
    return -1;
}

// hamming distance
inline double dist_hamming(const vector<int> &encode_1, const vector<int> &encode_2, const ReadRange &range_1, const ReadRange &range_2,
                           const vector<int> &var_cdf, vector<bool> &temp_array, int min_overlap = 500)
{
    ReadRange range_overlap(range_1.first >= range_2.first ? range_1.first : range_2.first, range_1.second <= range_2.second ? range_1.second : range_2.second);
    
    // get overlap length
    int overlap = range_overlap.second - range_overlap.first + 1;
    if (overlap < min_overlap)
        return -1;
    
    // get number of variants in overlap region
    int n_var = var_cdf[range_overlap.second] - var_cdf[range_overlap.first];
    if (n_var < 0)
        throw runtime_error("dist_hamming : n_var < 0");
    if (n_var == 0)
        return -2;
    
    // compare reads
    //int n_miss = (int)encode_1.size();
    
    int n_miss = 0;
    // read_1
    for (auto i = 0; i < encode_1.size(); ++i){
        if (encode_1[i] >= 4*range_overlap.first && encode_1[i] <= 4*range_overlap.second+3){
            temp_array[encode_1[i]] = true;
            ++n_miss;
        }
    }
    
    // read_2
    for (auto i = 0; i < encode_2.size(); ++i){
        if (encode_2[i] >= 4*range_overlap.first && encode_2[i] <= 4*range_overlap.second+3){
            if (temp_array[encode_2[i]])
                --n_miss;
            else
                ++n_miss;
        }
    }

    // clear up temp_array
    for (auto i = 0; i < encode_1.size(); ++i)
        temp_array[encode_1[i]] = false;
    
    return (double)n_miss / n_var;
}

inline void get_var_cdf(vector<int> &var_cdf, string var_file, size_t genome_size)
{
    ifstream fs_var_file;
    open_infile(fs_var_file, var_file);
    var_cdf.resize(genome_size, 0);
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
            var_cdf[i] = cur_count;
        
        ++cur_count;
        pre_pos = cur_pos + 1;
    }
    
    for (auto i = pre_pos; i < var_cdf.size(); ++i)
        var_cdf[i] = cur_count;
    
    fs_var_file.close();
}

inline size_t get_genome_size(const vector<ReadRange> &reads_range)
{
    size_t genome_size = 0;
    for (int i=0; i<(int)reads_range.size(); ++i)
        genome_size = genome_size < reads_range[i].second ? reads_range[i].second : genome_size;
    ++genome_size;

    return genome_size;
}

inline int get_nvar(const ReadRange &range_1, const ReadRange &range_2, const vector<int> &var_cdf)
{
    ReadRange range_overlap(range_1.first >= range_2.first ? range_1.first : range_2.first, range_1.second <= range_2.second ? range_1.second : range_2.second);
    
    // get overlap length
    int overlap = range_overlap.second - range_overlap.first + 1;
    if (overlap <= 0)
        return -1;
    
    // get number of variants in overlap region
    int n_var = var_cdf[range_overlap.second] - var_cdf[range_overlap.first];
    
    if (n_var < 0)
        throw runtime_error("dist_hamming : n_var < 0");
    
    if (n_var == 0)
        return -2;

    return n_var;
}
#endif
