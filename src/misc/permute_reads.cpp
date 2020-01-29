//
//  permute_reads.cpp
//  iGDA
//
//  Created by Zhixing Feng on 1/28/20.
//  Copyright © 2020 Zhixing Feng. All rights reserved.
//

#include "permute_reads.h"
#include "./io.h"

void permute_encodefile(string encode_file, string pileup_file, string outfile, int seed)
{
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    ifstream fs_pufile;
    
    // load pileup data
    unordered_map<int, vector<double> > pu_data;
    open_infile(fs_pufile, pileup_file);
    while (true) {
        string buf;
        getline(fs_pufile, buf);
        if (fs_pufile.eof())
            break;
        
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size() != 9)
            throw runtime_error("permute_encodefile(): buf_vec.size() != 9 in " + pileup_file);
        if (stod(buf_vec[8]) <= 0) continue;
        
        vector<double> cur_freq(4, -1);
        cur_freq[0] = stod(buf_vec[4]) / stod(buf_vec[8]);
        cur_freq[1] = stod(buf_vec[5]) / stod(buf_vec[8]);
        cur_freq[2] = stod(buf_vec[6]) / stod(buf_vec[8]);
        cur_freq[3] = stod(buf_vec[7]) / stod(buf_vec[8]);
        
        pu_data[stoi(buf_vec[0])] = cur_freq;
    }
    fs_pufile.close();
    
    // permute encode_data
    pcg32 rng(seed);
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for (auto i = 0; i < encode_data.size(); ++i){
        for (auto j = 0; j < encode_data[i].size(); ++j){
            // get substitution rate of the current locus
            int locus = int(encode_data[i][j] / 4);
            vector<double> cur_freq = pu_data[locus];
            double sub_rate = 0;
            for (auto k = 0; k < cur_freq.size(); ++k)
                sub_rate = cur_freq[k];
            if (sub_rate < 0 || sub_rate > 1)
                throw runtime_error("permute_encodefile(): sub_rate < 0 || sub_rate > 1");
            
            // toss a coin to determine if there is a substitution
            vector<double> sub_prob = {1 - sub_rate, sub_rate};
            vector<int> sub_cand = {0, 1};
            //cout << sub_prob << endl;
            //cout << sub_cand << endl;
            int is_sub = gen_binom(sub_prob, sub_cand, rng);
            
            if (is_sub == 1){
                vector<int> cand = {0,1,2,3};
                //cout << encode_data[i][j] << endl;
                int new_subtype = gen_binom(cur_freq, cand, rng);
                fs_outfile << 4*locus + new_subtype << "\t";
            }
        }
        fs_outfile << endl;
    }
    fs_outfile.close();
}
