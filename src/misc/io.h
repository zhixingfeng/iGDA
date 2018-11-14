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
#include "../modules/alignreader/alignreadersam.h"
#include <stxxl.h>
#include <zlib.h>
#include "../../tools/kseq.h"
#include "../../include/utils.h"
#include "./basic.h"

typedef stxxl::VECTOR_GENERATOR<vector<int64_t>, 16, 16, 10*1048576>::result stxxl_vector_type;
typedef stxxl::VECTOR_GENERATOR<vector<int>, 16, 16, 10*1048576>::result stxxl_vector_type_int;

typedef pair<int,int> ReadRange;

struct ReadMatch
{
    ReadMatch(): match_rate(0), n_overlap(0), start(0), end(0){}
    ReadMatch(vector<int> diff, vector<int> matches, double match_rate, int n_overlap, int read_id, int start, int end):
    diff(diff), matches(matches), match_rate(match_rate), n_overlap(n_overlap), read_id(read_id), start(start), end(end){}
    vector<int> diff;
    vector<int> matches;
    double match_rate;
    int n_overlap;
    int read_id;
    int start;
    int end;
};


// read cmpreads from .cmpreads files (text format)
inline bool loadcmpreads_txt(stxxl::vector<vector<int> > &cmpreads_data, string cmpreads_file)
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

// read cmpreads from .cmpreads files (binary format)
inline bool loadcmpreads(stxxl_vector_type_int &cmpreads_data, string cmpreads_file)
{
    // first scan cmpreads_file to get number of candidates
    FILE * p_cmpreads_file = fopen(cmpreads_file.c_str(), "rb");
    if (p_cmpreads_file == NULL)
        throw runtime_error("DForestSNVMax::run(): fail to open cmpreads_file");
    
    int64_t n_cand = 0;
    while(1){
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_cmpreads_file);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_cmpreads_file);
        if (feof(p_cmpreads_file))
            break;
        ++n_cand;
    }
    fclose(p_cmpreads_file);
    
    cout << "n_cand = " << n_cand << endl;
    // load the data
    cmpreads_data.resize(n_cand);
    //cmpreads_data = stxxl::vector<vector<int> >(n_cand, 1);
    
    p_cmpreads_file = fopen(cmpreads_file.c_str(), "rb");
    if (p_cmpreads_file == NULL)
        throw runtime_error("DForestSNVMax::run(): fail to open cmpreads_file");
    
    n_cand = 0;
    while(1){
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_cmpreads_file);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_cmpreads_file);
        if (feof(p_cmpreads_file))
            break;
        cmpreads_data[n_cand] = cand_loci;
        ++n_cand;
    }
    fclose(p_cmpreads_file);

    
    
    
    /*ofstream fs_outfile;
    open_outfile(fs_outfile, "./tmp_out_read.txt");
    while(1){
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_cmpreads_file);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_cmpreads_file);
        if (feof(p_cmpreads_file))
            break;
        fs_outfile << cand_loci << ',' << endl;
        cmpreads_data.push_back(cand_loci);
    }
    fs_outfile.close();
    fclose(p_cmpreads_file);
    
    open_outfile(fs_outfile, "./tmp_out_read_stxxl.txt");
    for (int64_t i = 0; i < cmpreads_data.size(); ++i)
        fs_outfile << cmpreads_data[i] << ',' << endl;
    fs_outfile.close();*/
    return true;
}


// read cmpreads_diff from .cmpreads_diff files (binary format)
inline bool loadcmpreadsdiff(stxxl::vector<ReadMatch> &cmpreadsdiff_data, string cmpreadsdiff_file)
{
    FILE *p_infile = fopen(cmpreadsdiff_file.c_str(), "rb");
    if (p_infile == NULL)
        runtime_error("fail to open cmpreadsdiff_file");
    
    while(true){
        // read line
        int read_id, start, end;
        int cand_loci_size, cand_loci_diff_size;
        
        fread(&read_id, sizeof(int), 1, p_infile);
        fread(&start, sizeof(int), 1, p_infile);
        fread(&end, sizeof(int), 1, p_infile);
        
        fread(&cand_loci_size, sizeof(int), 1, p_infile);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_infile);
        
        fread(&cand_loci_diff_size, sizeof(int), 1, p_infile);
        vector<int> cand_loci_diff(cand_loci_diff_size,-1);
        fread(&cand_loci_diff[0], sizeof(int), cand_loci_diff_size, p_infile);
        
        if (feof(p_infile))
            break;
        
        cmpreadsdiff_data.push_back(ReadMatch(cand_loci_diff, cand_loci, 0, 0, read_id, start, end));
        
    }
    fclose(p_infile);
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

// format: m = m5, s = sam, a = auto
inline bool loadreadsrange(vector<ReadRange> &reads_range, string align_file, char format='a')
{
    // setup alignreader
    AlignReader *p_alignreader;
    AlignReaderM5 alignreaderm5;
    AlignReaderSam alignreadersam;
    switch(format){
        case 'm':
            p_alignreader = &alignreaderm5;
            break;
        case 's':
            p_alignreader = &alignreadersam;
            break;
        case 'a':
            if (align_file.substr(align_file.size()-4,4) == ".sam")
                p_alignreader = &alignreadersam;
            else
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

inline void read_fasta(const string fasta_file, unordered_map<string, string> &fasta_data)
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

// get number of homopolymer blocks between each locus
inline vector<int64_t> get_homo_blocks(const string fasta_file)
{
    unordered_map<string, string> fasta_data;
    read_fasta(fasta_file, fasta_data);
    if (fasta_data.size()!=1)
        throw runtime_error("get_homo_blocks: fasta_data.size()!=1");
    auto it = fasta_data.begin();
    string cur_seq = it->second;
    if (cur_seq.size() == 0)
        throw runtime_error("get_homo_blocks: cur_seq.size() == 0");
    vector<int64_t> homo_blocks(cur_seq.size(), -1);
    char pre_base = '#';
    int64_t cur_idx = -1;
    for (auto i = 0; i < cur_seq.size(); ++i){
        char cur_base =  toupper(cur_seq[i]);
        if (cur_base != pre_base){
            pre_base = cur_base;
            ++cur_idx;
        }
        homo_blocks[i] = cur_idx;
    }
    
    return homo_blocks;
    
}


// load var_file
vector<int> load_varfile(string var_file)
{
    ifstream fs_var_file;
    open_infile(fs_var_file, var_file);

    vector<int> var_data;
    while (1){
        string buf;
        getline(fs_var_file, buf);
        if (fs_var_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in var_file");
        int cur_var = stod(buf_vec[2]);
        var_data.push_back(cur_var);
    }

    return var_data;
}


#endif /* io_h */
