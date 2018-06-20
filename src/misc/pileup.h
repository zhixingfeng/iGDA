//
//  pileup.h
//  iGDA
//
//  Created by Zhixing Feng on 16/10/3.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef pileup_h
#define pileup_h

#include "./cmpreads.h"
#include "basic.h"
#include <stxxl.h>
// coverage should NOT exceed range of int !!!!!
////////// pileup read ID and locus are 0-based ///////////
struct ConsensusSeq
{
    vector<int> pu_var_count;
    vector<int> pu_read_count;
    vector<double> prop;
    vector<int> cons_seq;
    vector<int> seed;
};

// pileup variants (input from memory/stxxl)
inline vector<vector<int> > pileup_var(const vector<vector<int> > &encode_data, int64_t &n_reads)
{
    // load encode data
    n_reads = encode_data.size();
    
    // calculate size of pileup vector
    int pu_size = 0;
    for (int i=0; i<(int)encode_data.size(); i++)
        for (int j=0; j<(int)encode_data[i].size(); j++)
            pu_size = encode_data[i][j] > pu_size ? encode_data[i][j] : pu_size;
    pu_size++;
    
    // pileup
    vector<vector<int> > pu(pu_size, vector<int>());
    for (int i=0; i<(int)encode_data.size(); i++)
        for (int j=0; j<(int)encode_data[i].size(); j++)
            pu[encode_data[i][j]].push_back(i);
    
    return pu;
}

// pileup variants (input from file)
inline vector<vector<int> > pileup_var(string encode_file, int64_t &n_reads)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    n_reads = encode_data.size();
    
    // calculate size of pileup vector
    int pu_size = 0;
    for (int i=0; i<(int)encode_data.size(); i++)
        for (int j=0; j<(int)encode_data[i].size(); j++)
            pu_size = encode_data[i][j] > pu_size ? encode_data[i][j] : pu_size;
    pu_size++;
    
    // pileup
    vector<vector<int> > pu(pu_size, vector<int>());
    for (int i=0; i<(int)encode_data.size(); i++)
        for (int j=0; j<(int)encode_data[i].size(); j++)
            pu[encode_data[i][j]].push_back(i);
    
    return pu;
}

inline int get_pu_var_size(const vector<vector<int> > &encode_data, const vector<int> &idx = vector<int>())
{
    int pu_size = 0;
    if (idx.size() > 0){
        for (int i : idx)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                pu_size = encode_data[i][j] > pu_size ? encode_data[i][j] : pu_size;
    }else{
        for (int i=0; i<(int)encode_data.size(); i++)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                pu_size = encode_data[i][j] > pu_size ? encode_data[i][j] : pu_size;
    }
    pu_size++;
    return pu_size;
}

// pileup variants (input from encode_data and reads_range)
inline vector<vector<int> > pileup_var(const vector<vector<int> > &encode_data, const vector<int> &idx = vector<int>())
{
    int pu_size = get_pu_var_size(encode_data, idx);
    
    // pileup
    vector<vector<int> > pu(pu_size, vector<int>());
    if (idx.size() > 0){
        for (int i : idx)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                pu[encode_data[i][j]].push_back(i);
    }else{
        for (int i=0; i<(int)encode_data.size(); i++)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                pu[encode_data[i][j]].push_back(i);
    }
    return pu;
}

inline void pileup_var_online(vector<vector<int> > &pu, const vector<int> &cur_encode_data, int cur_read_id)
{
    for (auto i = 0; i < cur_encode_data.size(); ++i)
        pu[cur_encode_data[i]].push_back(cur_read_id);
}

// pileup variants count
inline vector<int> pileup_var_count(const vector<vector<int> > &encode_data, const vector<int> &idx = vector<int>())
{
    int pu_size = get_pu_var_size(encode_data, idx);
    
    // pileup
    vector<int> pu_count(pu_size, 0);
    if (idx.size() > 0){
        for (int i : idx)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                ++pu_count[encode_data[i][j]];
    }else{
        for (int i=0; i<(int)encode_data.size(); i++)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                ++pu_count[encode_data[i][j]];
    }
    return pu_count;
}

// pileup count online
inline void pileup_var_online_count(vector<int> &pu_count, const vector<int> &cur_encode_data)
{
    for (auto i = 0; i < cur_encode_data.size(); ++i)
        ++pu_count[cur_encode_data[i]];
}

inline void pileup_var_online_count_pop(vector<int> &pu_count, const vector<int> &cur_encode_data)
{
    for (auto i = 0; i < cur_encode_data.size(); ++i)
        --pu_count[cur_encode_data[i]];
}


// pileup reads (input from memory/stxxl)
inline vector<vector<int> > pileup_reads_m5(const stxxl::vector<Align> &align_data, int64_t &n_reads)
{
    vector<ReadRange> reads_range;
    for (int i=0; i<(int)align_data.size(); ++i)
        reads_range.push_back(ReadRange(align_data[i].tStart, align_data[i].tEnd));
    
    n_reads = reads_range.size();
    
    // get size of pileup vector
    int pu_size=0;
    for (int i=0; i<(int)reads_range.size(); i++)
        pu_size = reads_range[i].second > pu_size ? reads_range[i].second : pu_size;
    pu_size++;
    
    // pileup
    vector<vector<int> > pu(pu_size, vector<int>());
    for (int i=0; i<(int)reads_range.size(); i++)
        for(int j=reads_range[i].first; j<=reads_range[i].second; j++)
            pu[j].push_back(i);
    return pu;
}

inline vector<vector<int> > pileup_reads(const stxxl::vector<Align> &align_data, int64_t &n_reads, char format = 'm')
{
    switch (format) {
        case 'm':
            return pileup_reads_m5(align_data, n_reads);
            break;
            
        default:
            return vector<vector<int> >(0,vector<int>(0,1));
            break;
    }
}

// pileup reads (input from file)
inline vector<vector<int> > pileup_reads_m5(string align_file, int64_t &n_reads, bool rm_del)
{
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file, 'm');
    n_reads = reads_range.size();
    
    // get size of pileup vector
    int pu_size=0;
    for (int i=0; i<(int)reads_range.size(); i++)
        pu_size = reads_range[i].second > pu_size ? reads_range[i].second : pu_size;
    pu_size++;
    
    // pileup
    vector<vector<int> > pu(pu_size, vector<int>());
    if (rm_del){
        // remove deletions
        AlignReaderM5 alignreaderm5;
        Align align;
        alignreaderm5.open(align_file);
        int k = 0;
        while(alignreaderm5.readline(align)){
            if (align.tAlignedSeq.size() != align.qAlignedSeq.size())
                throw runtime_error("pileup_reads_m5(): align.tAlignedSeq.size() != align.qAlignedSeq.size()");
            
            int cur_locus = align.tStart;
            for (int i = 0; i < (int)align.tAlignedSeq.size(); ++i){
                if (align.tAlignedSeq[i]!='-'){
                    if (align.qAlignedSeq[i]!='-')
                        pu[cur_locus].push_back(k);
                    ++cur_locus;
                }
            }
            
            ++k;
        }
        alignreaderm5.close();
        
    }else{
        // don't remove deletions
        for (int i=0; i<(int)reads_range.size(); i++)
            for(int j=reads_range[i].first; j<=reads_range[i].second; j++)
                pu[j].push_back(i);
    }
    return pu;
}

inline int get_pu_read_size(const vector<ReadRange> &reads_range, const vector<int> &idx = vector<int>())
{
    int pu_size=0;
    if (idx.size() > 0){
        for (int i : idx)
            pu_size = reads_range[i].second > pu_size ? reads_range[i].second : pu_size;
    }else{
        for (int i=0; i<(int)reads_range.size(); i++)
            pu_size = reads_range[i].second > pu_size ? reads_range[i].second : pu_size;
    }
    pu_size++;
    return pu_size;
}

inline vector<vector<int> > pileup_reads_m5(const vector<ReadRange> &reads_range, const vector<int> &idx = vector<int>())
{
    int pu_size = get_pu_read_size(reads_range, idx);
    vector<vector<int> > pu(pu_size, vector<int>());
    if (idx.size() > 0){
        for (int i : idx)
            for(int j=reads_range[i].first; j<=reads_range[i].second; j++)
                pu[j].push_back(i);
    }else{
        for (int i=0; i<(int)reads_range.size(); i++)
            for(int j=reads_range[i].first; j<=reads_range[i].second; j++)
                pu[j].push_back(i);
    }
    return pu;
}

inline void pileup_reads_m5_online(vector<vector<int> > &pu, const ReadRange &cur_reads_range, int cur_read_id)
{
    for(int i=cur_reads_range.first; i<=cur_reads_range.second; ++i)
        pu[i].push_back(cur_read_id);
}


// pileup reads count
inline vector<int> pileup_reads_m5_count(const vector<ReadRange> &reads_range, const vector<int> &idx = vector<int>())
{
    int pu_size = get_pu_read_size(reads_range, idx);
    vector<int> pu_count(pu_size, 0);
    if (idx.size() > 0){
        for (int i : idx)
            for(int j=reads_range[i].first; j<=reads_range[i].second; j++)
                ++pu_count[j];
    }else{
        for (int i=0; i<(int)reads_range.size(); i++)
            for(int j=reads_range[i].first; j<=reads_range[i].second; j++)
                ++pu_count[j];
    }
    return pu_count;
}

inline void pileup_reads_m5_online_count(vector<int> &pu_count, const ReadRange &cur_reads_range)
{
    for(int i=cur_reads_range.first; i<=cur_reads_range.second; ++i)
        ++pu_count[i];
}

inline void pileup_reads_m5_online_count_pop(vector<int> &pu_count, const ReadRange &cur_reads_range)
{
    for(int i=cur_reads_range.first; i<=cur_reads_range.second; ++i)
        --pu_count[i];
}


inline vector<vector<int> > pileup_reads(string align_file, int64_t &n_reads, bool rm_del = true, char format = 'm')
{
    switch (format) {
        case 'm':
            return pileup_reads_m5(align_file, n_reads, rm_del);
            break;
            
        default:
            return vector<vector<int> >(0,vector<int>(0,1));
            break;
    }
}

inline void print_pileup(const vector<vector<int> > &pileup_data, string outfile)
{
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for(int i=0; i<(int)pileup_data.size(); ++i)
        fs_outfile << pileup_data[i] << endl;
    fs_outfile.close();
}

inline vector<vector<int> > filter_pileup_var(const vector<vector<int> > &pu_var, const vector<vector<int> > &pu_read, int64_t n_reads)
{
    if (pu_var.size() > 4*pu_read.size()+3)
        throw runtime_error("filter_pileup_var(): pu_var.size() > 4*pu_read.size()");

    vector<vector<int> > pu_var_ft(pu_var.size(), vector<int>());
    vector<int> temp_vec(n_reads, -1);
    
    for (int i = 0; i < (int)pu_read.size(); ++i){
        if (i > int ((pu_var.size()-1) / 4))
            break;

        // fill in temp_vec
        for (int j = 0; j < (int)pu_read[i].size(); ++j)
            temp_vec[pu_read[i][j]] = i;
        
        // only keep variants covered by reads
        for (int k = 0; k <= 3; ++k){
            for (int j = 0; j < (int)pu_var[4*i+k].size(); ++j){
                if (temp_vec[pu_var[4*i+k][j]] == i)
                    pu_var_ft[4*i+k].push_back(pu_var[4*i+k][j]);
            }
        }
    }
    
    return pu_var_ft;
}

// get consensus
inline void get_consensus(ConsensusSeq &cons, const vector<int> &pu_var_count, const vector<int> &pu_read_count, int min_cvg = 20)
{
    if (floor(double(pu_var_count.size()-1) / 4) > pu_read_count.size() - 1)
        throw runtime_error("ann_clust: floor(double(pu_var_count.size()-1) / 4) > pu_read_count.size() - 1");

    cons.pu_var_count = pu_var_count;
    cons.pu_read_count = pu_read_count;
    cons.prop.resize(pu_var_count.size(),-1);
    for (auto i = 0; i < pu_var_count.size(); ++i){
        int i_r = int(i/4);
        if (pu_read_count[i_r] >= min_cvg)
              cons.prop[i] = (double) pu_var_count[i] / pu_read_count[i_r];
        else
              cons.prop[i] = -1;
        if (cons.prop[i] > 0.5)
            cons.cons_seq.push_back(i);
    }
}



#endif /* pileup_h */
