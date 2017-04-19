//
//  hclust.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/18.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "hclust.h"

void HClust::mask(string encode_file, string region_file, string out_file, bool is_0_based)
{
    if (ptr_aligncoder==NULL)
        throw runtime_error("ptr_aligncoder is NULL");
    
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, 0);

    // fill temp_array with encoded region
    ifstream fs_region;
    open_infile(fs_region, region_file);
    int cur_pos;
    char cur_base;
    while (1){
        fs_region >> cur_pos >> cur_base;
        if (fs_region.eof())
            break;
        if (!is_0_based){
            cur_pos = cur_pos - 1;
            if (cur_pos < 0)
                throw runtime_error("cur_pos < 0");
        }
        temp_array[ptr_aligncoder->binary_code(cur_pos, cur_base)] = 1;
    }
    fs_region.close();
    
    // scan encode data output matched part
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i = 0; i < encode_data.size(); i++){
        for (int j = 0; j < encode_data[i].size(); j++){
            if (temp_array[encode_data[i][j]]==1)
                fs_out << encode_data[i][j] << '\t';
        }
        fs_out << endl;
    }
    fs_out.close();
    
}

void HClust::dist(string encode_file, string align_file, string out_file, bool is_nmiss)
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
            int n_overlap = (reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second) -
            (reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first) + 1;
            
            if (n_overlap <= 0){
                fs_out << -1 << '\t';
                continue;
            }
            int n_miss = (int)encode_data[i].size();
            for (int k = 0; k < encode_data[j].size(); k++){
                if (temp_array[encode_data[j][k]] == i)
                    n_miss--;
                else
                    n_miss++;
            }
            if (n_miss < 0)
                throw runtime_error("n_miss < 0");
            if (is_nmiss)
                fs_out << n_miss << '\t';
            else
                fs_out << (double)n_miss / n_overlap << '\t';
        }
        fs_out << endl;
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();
}







