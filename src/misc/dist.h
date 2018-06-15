//
//  dist.h
//  iGDA
//
//  Created by Zhixing Feng on 6/15/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#ifndef iGDA_dist_h
#define iGDA_dist_h

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
    int n_miss = (int)encode_1.size();
    
    // read_1
    for (auto i = 0; i < encode_1.size(); ++i)
        temp_array[encode_1[i]] = true;
    
    // read_2
    for (auto i = 0; i < encode_2.size(); ++i){
        if (temp_array[encode_2[i]])
            --n_miss;
        else
            ++n_miss;
    }
    
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
        ++cur_count;
        for (int i=pre_pos; i<=cur_pos; ++i)
            var_cdf[i] = cur_count;
        pre_pos = cur_pos + 1;
    }
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


#endif