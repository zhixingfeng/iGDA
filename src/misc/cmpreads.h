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
    cout << "load encode" << endl;
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    cout << "load reads range" << endl;
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
        // hash the reads to be compared
        G_Hash cur_variant;
        for (int j=0; j<(int)encode_data[i].size(); j++)
            cur_variant.insert(pair<int,int>(encode_data[i][j],1));
        
        // compare other reads to cur_variant
        for (int j=i+1; j<(int)encode_data.size(); j++){
            p_out_file << "code:";
            for (int k=0; k<(int)encode_data[j].size();k++){
                if (cur_variant.find(encode_data[j][k]) != cur_variant.end()){
                    p_out_file << encode_data[j][k] << ',';
                }
            }
            p_out_file << '\t' << reads_range[i].first << ',' << reads_range[i].second << '\t';
            p_out_file << reads_range[j].first << ',' << reads_range[j].second << endl;
        }
    }
    cout << encode_data.size() << endl;
    p_out_file.close();
    
    return true;
}



#endif /* cmpreads_h */
