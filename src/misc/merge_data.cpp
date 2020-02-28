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
