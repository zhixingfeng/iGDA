//
//  cmpreads.h
//  iGDA
//
//  Created by Zhixing Feng on 16/9/7.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.

//  pairwise compare encoded reads to find consistently occuring variants

#ifndef cmpreads_h
#define cmpreads_h
#include "../../include/headers.h"
#include "../modules/alignreader/alignreaderm5.h"
#include "../../tools/khash.h"

KHASH_MAP_INIT_INT(32, char)
//typedef map<int,int> G_Hash;
typedef unordered_map<int,int> G_Hash;
typedef pair<int,int> ReadRange;

inline bool loadencodedata(vector<vector<int> > &encode_data, string encode_file)
{
    ifstream p_encode_file; open_infile(p_encode_file, encode_file);
    while (true) {
        string buf;
        getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        encode_data.push_back(split_int(buf, '\t'));
    }
    p_encode_file.close();
    return true;
}

// format: m=m5
inline bool loadreadsrange(vector<ReadRange> &reads_range, string align_file, char format='m')
{
    // setup alignreader
    AlignReader *p_alignreader;
    AlignReaderM5 alignreaderm5;
    switch(format){
        case 'm':
            p_alignreader = &alignreaderm5;
            break;
        default:
            throw runtime_error("loadreadsranges: unsupported format.");
    }
    
    // load alignment data
    Align align;
    p_alignreader->open(align_file);
    while(p_alignreader->readline(align))
        reads_range.push_back(ReadRange(align.tStart, align.tEnd));
    
    p_alignreader->close();
    
    return true;
}

inline bool cmpreads(string encode_file, string align_file, string out_file)
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
    
    // pairwise comparison
    ofstream p_out_file;
    open_outfile(p_out_file, out_file);
    for (int i=0; i<(int)(encode_data.size()-1); i++){
        if ((i+1)%100==0) cout << i+1 << '\r';
        // khash the reads to be compared
        khiter_t k_it; int ret;
        khash_t(32) *cur_variant = kh_init(32);
        for (int j=0; j<(int)encode_data[i].size(); j++)
            k_it = kh_put(32, cur_variant, encode_data[i][j], &ret);
        
        // compare other reads to cur_variant
        int n_max = 0;
        int idx_max = 0;
        vector<int> max_match;
        for (int j=0; j<(int)encode_data.size(); j++){
            if (j==i) continue;
            if (reads_range[i].first >= reads_range[j].second || reads_range[j].first >= reads_range[i].second)
                continue;
            
            vector<int> cur_match;
            for (int k=0; k<(int)encode_data[j].size(); k++){
                if (kh_get(32, cur_variant, encode_data[j][k]) != kh_end(cur_variant))
                    cur_match.push_back(encode_data[j][k]);
            }
            if (cur_match.size() > n_max){
                n_max = (int)cur_match.size();
                idx_max = j;
                max_match = cur_match;
            }
            
            //p_out_file << i+1 << ',' << j+1 << "\t";
            //for (int k=0; k<(int)cur_match.size(); k++)
            //    p_out_file << cur_match[k] << ',';
            //p_out_file << '\t' << reads_range[i].first << ',' << reads_range[i].second << '\t';
            //p_out_file << reads_range[j].first << ',' << reads_range[j].second << endl;
        }
        kh_destroy(32, cur_variant);
        
        // output max_match for read i
        if (n_max==0) continue;
        p_out_file << i+1 << ',' << idx_max+1 << "\t";
        for (int k=0; k<(int)max_match.size(); k++)
            p_out_file << max_match[k] << ',';
        p_out_file << '\t' << reads_range[i].first << ',' << reads_range[i].second << '\t';
        p_out_file << reads_range[idx_max].first << ',' << reads_range[idx_max].second << endl;
    }
    cout << encode_data.size() << endl;
    p_out_file.close();
    
    return true;
}

/*inline bool cmpreads(string encode_file, string align_file, string out_file)
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
    
    // pairwise comparison
    ofstream p_out_file;
    open_outfile(p_out_file, out_file);
    for (int i=0; i<(int)(encode_data.size()-1); i++){
        if ((i+1)%100==0) cout << i+1 << '\r';
        // khash the reads to be compared
        khiter_t k_it; int ret;
        khash_t(32) *cur_variant = kh_init(32);
        for (int j=0; j<(int)encode_data[i].size(); j++)
            k_it = kh_put(32, cur_variant, encode_data[i][j], &ret);
        
        // compare other reads to cur_variant
        for (int j=i+1; j<(int)encode_data.size(); j++){
            if (reads_range[i].first >= reads_range[j].second || reads_range[j].first >= reads_range[i].second)
                continue;
            
            
            vector<int> cur_match;
            for (int k=0; k<(int)encode_data[j].size();k++){
                if (kh_get(32, cur_variant, encode_data[j][k]) != kh_end(cur_variant))
                    cur_match.push_back(encode_data[j][k]);
            }
            if (cur_match.size()==0)
                continue;
            p_out_file << i+1 << ',' << j+1 << "\t";
            for (int k=0; k<(int)cur_match.size(); k++)
                p_out_file << cur_match[k] << ',';
            p_out_file << '\t' << reads_range[i].first << ',' << reads_range[i].second << '\t';
            p_out_file << reads_range[j].first << ',' << reads_range[j].second << endl;
        }
        
        kh_destroy(32, cur_variant);
        
    }
    cout << encode_data.size() << endl;
    p_out_file.close();
    
    return true;
}
*/


#endif /* cmpreads_h */
