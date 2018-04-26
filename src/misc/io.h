//
//  io.h
//  iGDA
//
//  Created by Zhixing Feng on 16/10/9.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef io_h
#define io_h
#include "../../include/headers.h"
#include "../modules/alignreader/alignreaderm5.h"
#include <stxxl.h>
#include <zlib.h>
#include "../../tools/kseq.h"
#include "../../include/utils.h"

typedef pair<int,int> ReadRange;

// read cmpreads from .cmpreads files (text format)
inline bool loadcmpreads(stxxl::vector<vector<int> > &cmpreads_data, string cmpreads_file)
{
    ifstream p_cmpreads_file; open_infile(p_cmpreads_file, cmpreads_file);
    while (true) {
        string buf;
        getline(p_cmpreads_file, buf);
        if (p_cmpreads_file.eof())
            break;
        cmpreads_data.push_back(split_int(buf, ','));
    }
    p_cmpreads_file.close();
    return true;

}

// load encode data
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

// idx must be increasingly ordered
inline void select_lines(const vector<int> &idx, string infile, string outfile)
{
    ifstream fs_infile; ofstream fs_outfile;
    open_infile(fs_infile, infile);
    open_outfile(fs_outfile, outfile);
    int i = 0;
    int k = 0;
    while(true){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof())
            break;
        if (k == idx[i]){
            fs_outfile << buf << endl;
            ++i;
        }
        if (i >= idx.size())
            break;
        ++k;
    }
    fs_infile.close();
    fs_outfile.close();
}

KSEQ_INIT(gzFile, gzread)

inline void read_fasta(string fasta_file, unordered_map<string, string> &fasta_data)
{
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fasta_file.c_str(), "r");
    if (fp == NULL)
        throw runtime_error("fail to open " + fasta_file);
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        auto it = fasta_data.find(seq->name.s);
        if (it != fasta_data.end())
            throw runtime_error("duplicated seqence " + string(seq->name.s));
        fasta_data[string(seq->name.s)] = string(seq->seq.s);
       
    }
    
    kseq_destroy(seq);
    gzclose(fp);
}

inline unordered_map<int, int> load_ncread(string check_follower_file, int min_follower = 5)
{
    unordered_map<int, int> ncread_id;
    
    ifstream fs_infile;
    open_infile(fs_infile, check_follower_file);
    while(true){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof())
            break;
        
        vector<string> buf_vec = split(buf, '\t');
        
        if (buf_vec.size() < 2)
            throw runtime_error("load_ncread(): buf_vec.size() < 2");
        
        int read_id = stoi(buf_vec[0]);
        int n_follower = stoi(buf_vec[1]);
        
        if (n_follower >= min_follower)
            ncread_id[read_id] = n_follower;
    }
    fs_infile.close();
    return ncread_id;
}



#endif /* io_h */
