//
//  cmpreads.h
//  iGDA
//
//  Created by Zhixing Feng on 16/9/7.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.

//  pairwise compare encoded reads to find consistently occuring variants

#ifndef cmpreads_h
#define cmpreads_h

#include "io.h"

inline bool cmpreads(string encode_file, string align_file, string out_file, double min_overlap = 0.25,
                     bool is_rm_single=true, bool is_binary=true)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (encode_data.size() != reads_range.size())
        throw runtime_error("cmpreads: size of encode_data and reads_range do not match.");
    
    cout << encode_data.size() << endl;
    
    // get the right-most variant location to determing size of template array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    FILE *p_out_binfile = NULL;
    ofstream p_out_textfile;
    if (is_binary)
        p_out_binfile = fopen(out_file.c_str(), "wb");
    else
        open_outfile(p_out_textfile, out_file);
    
    for (int i=0; i<(int)(encode_data.size()-1); i++){
        if ((i+1)%100==0) cout << i+1 << '\r';
        
        // fill the template array by the variants in the ith read
        for (int j = 0; j < encode_data[i].size(); j++)
            temp_array[encode_data[i][j]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=i+1; j<(int)encode_data.size(); j++){
            //if (reads_range[i].first > reads_range[j].second || reads_range[j].first > reads_range[i].second)
            //    continue;
            
            int n_overlap = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second -
                            reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first + 1;
            if (n_overlap < min_overlap * (reads_range[i].second - reads_range[i].first + 1) && n_overlap < min_overlap * (reads_range[j].second - reads_range[j].first + 1))
                continue;
            vector<int> cur_match;
            for (int k = 0; k < encode_data[j].size(); k++)
                if (temp_array[encode_data[j][k]] == i)
                    cur_match.push_back(encode_data[j][k]);
            
            if (is_rm_single){
                if (cur_match.size() < 2)
                    continue;
            }else{
                if (cur_match.size() == 0)
                    continue;
            }
            if (is_binary){
                int cur_match_size = (int)cur_match.size();
                fwrite(&cur_match_size, sizeof(int), 1, p_out_binfile);
                fwrite(&cur_match[0], sizeof(int), cur_match_size, p_out_binfile);
            }else{
                p_out_textfile << i+1 << ',' << j+1 << "\t";
                for (int k=0; k<(int)cur_match.size(); k++)
                    p_out_textfile << cur_match[k] << ',';
                p_out_textfile << '\t' << reads_range[i].first << ',' << reads_range[i].second << '\t';
                p_out_textfile << reads_range[j].first << ',' << reads_range[j].second << endl;
            }
        }
                
    }
    cout << encode_data.size() << endl;
    
    if (is_binary)
        fclose(p_out_binfile);
    else
        p_out_textfile.close();
    
    return true;
}



#endif /* cmpreads_h */
