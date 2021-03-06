//
//  merge_data.cpp
//  iGDA
//
//  Created by Zhixing Feng on 2020/2/27.
//  Copyright © 2020年 Zhixing Feng. All rights reserved.
//

#include "./merge_data.h"
void merge_encode(string m5_fofn_file, string encode_fofn_file, string out_encode_file)
{
    /*------------- load fofn files --------------*/
    // load m5 files
    vector<string> m5_files;
    ifstream p_m5_fofn_file; open_infile(p_m5_fofn_file, m5_fofn_file);
    while (true) {
        string buf;
        getline(p_m5_fofn_file, buf);
        if (p_m5_fofn_file.eof()) break;
        m5_files.push_back(buf);
    }
    p_m5_fofn_file.close();
    
    // load encode files
    vector<string> encode_files;
    ifstream p_encode_fofn_file; open_infile(p_encode_fofn_file, encode_fofn_file);
    while (true) {
        string buf;
        getline(p_encode_fofn_file, buf);
        if (p_encode_fofn_file.eof()) break;
        encode_files.push_back(buf);
    }
    p_encode_fofn_file.close();
    
    if (m5_files.size() != encode_files.size())
    throw runtime_error("merge_encode(): m5_files.size() != encode_files.size()");
    
    /*----------- merge encode file ------------*/
    unordered_map<string, set<int> > map_encode;
    for (auto i = 0; i < m5_files.size(); ++i){
        cout << "m5_file = " << m5_files[i] << endl;
        cout << "encode_file = " << encode_files[i] << endl;
        // load m5 data
        unordered_map<string, ReadRange> m5_data;
        vector<string> read_names;
        loadm5data(m5_data, read_names, m5_files[i]);
        
        // load encode data
        vector<vector<int> > encode_data;
        loadencodedata(encode_data, encode_files[i]);
        
        // check data
        if (m5_data.size() != encode_data.size() || read_names.size() != encode_data.size())
            throw runtime_error("merge_encode(): reads_range.size() != encode_data.size()");
        
        // merge
        for (auto j = 0; j < read_names.size() ; ++j){
            if (map_encode.find(read_names[j]) == map_encode.end())
                map_encode[read_names[j]] = set<int>();
            
            for (auto k = 0; k < encode_data[j].size(); ++k){
                if (encode_data[j][k] / 4 < m5_data[read_names[j]].first || encode_data[j][k] / 4 > m5_data[read_names[j]].second)
                    throw runtime_error("merge_encode(): unmatched encode_data and m5_data");
                map_encode[read_names[j]].insert(encode_data[j][k]);
            }
        }
    }
    
    /*------------- output -------------*/
    ofstream p_outfile; ofstream p_outfile_name;
    open_outfile(p_outfile, out_encode_file);
    open_outfile(p_outfile_name, out_encode_file + ".readname");
    for (auto it_i = map_encode.begin(); it_i != map_encode.end(); ++it_i){
        for (auto it_j = it_i->second.begin(); it_j != it_i->second.end(); ++it_j){
            p_outfile << *it_j << '\t';
        }
        p_outfile << endl;
        p_outfile_name << it_i->first << endl;
    }
    p_outfile.close();
    p_outfile_name.close();
}

void merge_m5(string m5_fofn_file, string readname_file, string out_m5_file)
{
    // load m5 files
    vector<string> m5_files;
    ifstream p_m5_fofn_file; open_infile(p_m5_fofn_file, m5_fofn_file);
    while (true) {
        string buf;
        getline(p_m5_fofn_file, buf);
        if (p_m5_fofn_file.eof()) break;
        m5_files.push_back(buf);
    }
    p_m5_fofn_file.close();
    
    // load readname
    vector<string> readnames;
    ifstream p_readname_file; open_infile(p_readname_file, readname_file);
    while (true) {
        string buf;
        getline(p_readname_file, buf);
        if (p_readname_file.eof()) break;
        readnames.push_back(buf);
    }
    p_readname_file.close();
    
    
    // merge m5 files
    unordered_map<string, ReadRange> m5_data;
    for (auto i = 0; i < m5_files.size(); ++i){
        cout << m5_files[i] << endl;
        unordered_map<string, ReadRange> cur_m5_data;
        vector<string> read_names;
        loadm5data(cur_m5_data, read_names, m5_files[i]);
        for (auto it = cur_m5_data.begin(); it != cur_m5_data.end(); ++it){
            auto it_hit = m5_data.find(it->first);
            if (it_hit == m5_data.end()){
                m5_data[it->first] = it->second;
            }else{
                if (it->second.first < it_hit->second.first)
                    it_hit->second.first = it->second.first;
                if (it->second.second > it_hit->second.second)
                    it_hit->second.second = it->second.second;
            }
        }
    }
    
    if (m5_data.size() != readnames.size())
        throw runtime_error("merge_m5(): m5_data.size() != readnames.size()");
    
    // print merged m5
    ofstream p_out_m5_file; open_outfile(p_out_m5_file, out_m5_file);
    for (auto i = 0; i < readnames.size(); ++i){
        auto it = m5_data.find(readnames[i]);
        if (it == m5_data.end())
            throw runtime_error("readname and m5_data do not match.");
        ReadRange cur_readrange = it->second;
        p_out_m5_file << it->first << ' ';
        for (auto k = 1; k <= 6; ++k)
            p_out_m5_file << "0 ";
        p_out_m5_file << cur_readrange.first << ' ' << cur_readrange.second + 1 << ' ';
        for (auto k = 1; k <= 9; ++k)
            p_out_m5_file << "0 ";
        p_out_m5_file<< "0" << endl;
    }
    p_out_m5_file.close();
}
