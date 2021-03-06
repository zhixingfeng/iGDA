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
#include "../modules/alignreader/alignreadersam.h"
//#include <stxxl.h>
// coverage should NOT exceed range of int !!!!!
////////// pileup read ID and locus are 0-based ///////////
struct ConsensusSeq
{
    ConsensusSeq(){rr_null = -1; rr_ind = -1; log_bf_null = -1000;  log_bf_ind = -1000; contig_count = 0; contig_cvg = 0;}
    vector<int> pu_var_count;
    vector<int> pu_read_count;
    vector<double> prop;
    vector<int> cons_seq;
    vector<int> seed;
    vector<int> neighbors_id;
    
    vector<int> tested_loci; // store test loci which pass the minimal coverage
    
    double contig_count;
    double contig_cvg;
    vector<int64_t> nn_reads_id;
    
    int start;
    int end;

    double log_bf_null; // test against the consensus
    double log_bf_ind; // test if variants are independent
    
    double rr_null; // relative risk against the consensus
    double rr_ind; // relative risk to test if variants are independent

};

inline void print_pileup_qv_count(const vector<vector<double> > &pu_qv_count, string outfile)
{
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for (auto i = 0; i < pu_qv_count.size(); ++i){
        // print locus and code
        fs_outfile << i << '\t' << int64_t(i/4) << '\t';
        // print read base
        char base = 'N';
        switch(i%4){
            case 0:
                base = 'A';
                break;
            case 1:
                base = 'C';
                break;
            case 2:
                base = 'G';
                break;
            case 3:
                base = 'T';
                break;
            default:
                break;
        }
        // print ref base
        char ref;
        switch(int(pu_qv_count[i][4])){
            case 1:
                ref = 'A';
                break;
            case 2:
                ref = 'C';
                break;
            case 3:
                ref = 'G';
                break;
            case 4:
                ref = 'T';
                break;
            case 5:
                ref = 'N';
                break;
            default:
                ref = '*';
                //throw runtime_error("base show be A, C, G, T or N");
                break;
        }
        
        // print base, ref and qv
        fs_outfile << base << '\t' << ref << '\t' << pu_qv_count[i][0] << '\t' << pu_qv_count[i][1] << '\t' << pu_qv_count[i][2] << '\t' << pu_qv_count[i][3] ;
        
        fs_outfile << endl;
    }
    
    fs_outfile.close();
}

inline void print_pileup_qv(const vector<vector<pair<int64_t, double> > > &pu_qv, string outfile)
{
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for (auto i = 0; i < pu_qv.size(); ++i){
        fs_outfile << i << '\t' << int64_t(i/4) << '\t';
        char base = 'N';
        switch(i%4){
            case 0:
                base = 'A';
                break;
            case 1:
                base = 'C';
                break;
            case 2:
                base = 'G';
                break;
            case 3:
                base = 'T';
                break;
            default:
                break;
        }
        fs_outfile << base << '\t';
        
        if (pu_qv[i].size() == 0){
            fs_outfile << -1 << '\t' << -1;
        }else{
            for (auto j = 0; j < pu_qv[i].size(); ++j)
                fs_outfile << pu_qv[i][j].second << ',';
            
            fs_outfile << '\t';
            
            for (auto j = 0; j < pu_qv[i].size(); ++j)
                fs_outfile << pu_qv[i][j].first << ',';
            
        }
        fs_outfile << endl;
    }
    
    fs_outfile.close();
}

inline vector<vector<double> > pileup_qv_count(const string sam_file, const string ref_fafile, bool is_var = true)
{
    AlignReaderSam alignreader;
    
    // get reference genome
    alignreader.getref(ref_fafile);
    
    // first scan the sam/bam file to get genome_size
    alignreader.open(sam_file);
    size_t g_size = 0;
    Align align;
    while(alignreader.readline(align)){
        if (align.tEnd + 1 > g_size)
            g_size = align.tEnd + 1;
    }
    alignreader.close();
    
    // initialize pileup
    
    vector<vector<double> > pu_qv_count(4*g_size, vector<double>(5,0)); // [0]="total qv", [1]="total var", [2]="mean qv", [3]="coverage(no indel)", [4]="reference, 1=A, 2=C, 3=G, 4=T, 5=N"
    
    // scan the sam/bam file again to pileup
    alignreader.open(sam_file);
    int64_t read_id = 0;
    while(alignreader.readline(align)){
        for (auto i = 0; i < align.qv.size(); ++i){
            // skip indels
            if (align.qAlignedSeq[i] == '-' || align.tAlignedSeq[i] == '-')
                continue;
            
            // calculate coverage
            ++pu_qv_count[4*align.qv_locus[i]][3];
            ++pu_qv_count[4*align.qv_locus[i] + 1][3];
            ++pu_qv_count[4*align.qv_locus[i] + 2][3];
            ++pu_qv_count[4*align.qv_locus[i] + 3][3];
            
            // get reference seq
            switch(align.tAlignedSeq[i]){
                case 'A':
                    for (auto k = 0; k <= 3; ++k)
                        pu_qv_count[4*align.qv_locus[i] + k][4] = 1;
                    break;
                case 'C':
                    for (auto k = 0; k <= 3; ++k)
                        pu_qv_count[4*align.qv_locus[i] + k][4] = 2;
                    break;
                case 'G':
                    for (auto k = 0; k <= 3; ++k)
                        pu_qv_count[4*align.qv_locus[i] + k][4] = 3;
                    break;
                case 'T':
                    for (auto k = 0; k <= 3; ++k)
                        pu_qv_count[4*align.qv_locus[i] + k][4] = 4;
                    break;
                case 'N':
                    for (auto k = 0; k <= 3; ++k)
                        pu_qv_count[4*align.qv_locus[i] + k][4] = 5;
                    break;
                default:
                    for (auto k = 0; k <= 3; ++k)
                        pu_qv_count[4*align.qv_locus[i] + k][4] = -1;
                    //throw runtime_error("base show be A, C, G, T or N");
                    break;
            }
            
            // pileup qv
            if (is_var && align.qAlignedSeq[i] == align.tAlignedSeq[i])
                continue;
            
            switch(align.qAlignedSeq[i]){
                case 'A':
                    pu_qv_count[4*align.qv_locus[i]][0] +=align.qv[i];
                    ++pu_qv_count[4*align.qv_locus[i]][1];
                    break;
                case 'C':
                    pu_qv_count[4*align.qv_locus[i] + 1][0] +=align.qv[i];
                    ++pu_qv_count[4*align.qv_locus[i] + 1][1];
                    break;
                case 'G':
                    pu_qv_count[4*align.qv_locus[i] + 2][0] +=align.qv[i];
                    ++pu_qv_count[4*align.qv_locus[i] + 2][1];
                    break;
                case 'T':
                    pu_qv_count[4*align.qv_locus[i] + 3][0] +=align.qv[i];
                    ++pu_qv_count[4*align.qv_locus[i] + 3][1];
                    break;
                case 'N':
                    break;
                default:
                    throw runtime_error("base show be A, C, G, T or N");
                    break;
            }
            
        }
        ++read_id;
    }
    alignreader.close();
    
    for (auto i = 0; i < pu_qv_count.size(); ++i)
        if (pu_qv_count[i][1] > 0)
            pu_qv_count[i][2] = (double)pu_qv_count[i][0] / pu_qv_count[i][1];
        else
            pu_qv_count[i][2] = 0;
    
    return pu_qv_count;
}

inline vector<vector<pair<int64_t, double> > > pileup_qv(const string sam_file, const string ref_fafile, bool is_var = true)
{
    AlignReaderSam alignreader;
    
    // get reference genome
    alignreader.getref(ref_fafile);
    
    // first scan the sam/bam file to get genome_size
    alignreader.open(sam_file);
    size_t g_size = 0;
    Align align;
    while(alignreader.readline(align)){
        if (align.tEnd + 1 > g_size)
            g_size = align.tEnd + 1;
    }
    alignreader.close();
    
    // initialize pileup
    
    vector<vector<pair<int64_t, double> > > pu_qv(4*g_size, vector<pair<int64_t, double> >());
    
    // scan the sam/bam file again to pileup
    alignreader.open(sam_file);
    int64_t read_id = 0;
    while(alignreader.readline(align)){
        for (auto i = 0; i < align.qv.size(); ++i){
            // skip indels
            if (align.qAlignedSeq[i] == '-' || align.tAlignedSeq[i] == '-')
                continue;
            
            // pileup
            if (is_var && align.qAlignedSeq[i] == align.tAlignedSeq[i])
                continue;
            
            switch(align.qAlignedSeq[i]){
                case 'A':
                    pu_qv[4*align.qv_locus[i]].push_back( pair<int64_t, double>(read_id, align.qv[i]) );
                    break;
                case 'C':
                    pu_qv[4*align.qv_locus[i] + 1].push_back( pair<int64_t, double>(read_id, align.qv[i]) );
                    break;
                case 'G':
                    pu_qv[4*align.qv_locus[i] + 2].push_back( pair<int64_t, double>(read_id, align.qv[i]) );
                    break;
                case 'T':
                    pu_qv[4*align.qv_locus[i] + 3].push_back( pair<int64_t, double>(read_id, align.qv[i]) );
                    break;
                case 'N':
                    break;
                default:
                    throw runtime_error("base show be A, C, G, T or N");
                    break;
            }
            
        }
        ++read_id;
    }
    alignreader.close();
    
    return pu_qv;
}

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
    int pu_size = -1;
    if (idx.size() > 0){
        for (int i : idx)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                pu_size = encode_data[i][j] > pu_size ? encode_data[i][j] : pu_size;
    }else{
        for (int i=0; i<(int)encode_data.size(); i++)
            for (int j=0; j<(int)encode_data[i].size(); j++)
                pu_size = encode_data[i][j] > pu_size ? encode_data[i][j] : pu_size;
    }
    //pu_size++;
    if (pu_size == -1) return 0;
    return 4*(pu_size/4) + 3 + 1;
}

// pileup variants (input from encode_data and reads_range)
inline vector<vector<int> > pileup_var(const vector<vector<int> > &encode_data, const vector<int> &idx = vector<int>())
{
    int pu_size = get_pu_var_size(encode_data, idx) + 4;
    
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
inline void pileup_var_online_count(vector<int> &pu_count, const vector<int> &cur_encode_data, unordered_set<int64_t> &mod_idx )
{
    for (auto i = 0; i < cur_encode_data.size(); ++i){
        ++pu_count[cur_encode_data[i]];
        mod_idx.insert(cur_encode_data[i]);
    }
}

inline void pileup_var_online_count_pop(vector<int> &pu_count, const vector<int> &cur_encode_data)
{
    for (auto i = 0; i < cur_encode_data.size(); ++i){
        --pu_count[cur_encode_data[i]];
        if (pu_count[cur_encode_data[i]] < 0)
            throw runtime_error("pu_count[cur_encode_data[i]] < 0");
    }
}


// pileup reads (input from memory/stxxl)
inline vector<vector<int> > pileup_reads_m5(const vector<Align> &align_data, int64_t &n_reads)
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

inline vector<vector<int> > pileup_reads(const vector<Align> &align_data, int64_t &n_reads, char format = 'a')
{
    switch (format) {
        case 'a':
            return pileup_reads_m5(align_data, n_reads);
            break;
            
        default:
            return vector<vector<int> >(0,vector<int>(0,1));
            break;
    }
}

// pileup reads only consider detected SNVs
inline vector<vector<int> > pileup_reads_m5_reduce(string align_file, string var_file, int64_t &n_reads)
{
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    n_reads = reads_range.size();
    
    // load detected var (coded)
    vector<int> detected_var = load_varfile(var_file);
    
    // get size of pileup vector
    int pu_size = -1;
    for (int i=0; i<(int)reads_range.size(); i++)
        pu_size = reads_range[i].second > pu_size ? reads_range[i].second : pu_size;
    pu_size++;
    
    if (pu_size <= 0) throw runtime_error("pileup_reads_m5_reduce: pu_size <= 0");
    
    // check detected locus
    vector<bool> is_detected(pu_size, false);
    for (auto i = 0; i < detected_var.size(); ++i){
        int cur_locus = int(detected_var[i] / 4);
        if (cur_locus >= pu_size) throw runtime_error("pileup_reads_m5_reduce: cur_locus >= pu_size");
        is_detected[cur_locus] = true;
    }
    
    // pileup
    vector<vector<int> > pu(pu_size, vector<int>());
    for (int i=0; i<(int)reads_range.size(); i++){
        for(int j=reads_range[i].first; j<=reads_range[i].second; j++){
            if (is_detected[j]){
                pu[j].push_back(i);
            }
        }
    }
    return pu;
}

// pileup reads (input from file)
inline vector<vector<int> > pileup_reads_m5(string align_file, int64_t &n_reads, bool rm_del = false)
{
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    n_reads = reads_range.size();
    
    // get size of pileup vector
    int pu_size = -1;
    for (int i=0; i<(int)reads_range.size(); i++)
        pu_size = reads_range[i].second > pu_size ? reads_range[i].second : pu_size;
    pu_size++;
    if (pu_size <= 0) throw runtime_error("pileup_reads_m5_reduce: pu_size <= 0");
    
    // pileup
    vector<vector<int> > pu(pu_size, vector<int>());
    if (rm_del){
        cout << "rm_del" << endl;
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
    int pu_size=-1;
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

inline void pileup_reads_m5_online_count(vector<int> &pu_count, const ReadRange &cur_reads_range, ReadRange &mod_range)
{
    for(int i=cur_reads_range.first; i<=cur_reads_range.second; ++i){
        ++pu_count[i];
        if (cur_reads_range.first < mod_range.first)
            mod_range.first = cur_reads_range.first;
        if (cur_reads_range.second > mod_range.second)
            mod_range.second = cur_reads_range.second;
    }
}

inline void pileup_reads_m5_online_count_pop(vector<int> &pu_count, const ReadRange &cur_reads_range)
{
    for(int i=cur_reads_range.first; i<=cur_reads_range.second; ++i){
        --pu_count[i];
        if (pu_count[i] < 0)
            throw runtime_error("pu_count[i] < 0");
    }
}


inline vector<vector<int> > pileup_reads(string align_file, int64_t &n_reads, char format = 'a', bool rm_del = false)
{
    if (rm_del)
        cout << "rm_del" << endl;
    switch (format) {
        case 'a':
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
    if (pu_var.size()-1 > 4*(pu_read.size()-1)+3){
        cout << "pu_var.size() = " << pu_var.size() << endl;
        cout << "pu_read.size() = " << pu_read.size() << endl;
        throw runtime_error("filter_pileup_var(): pu_var.size() > 4*pu_read.size()");
    }

    vector<vector<int> > pu_var_ft(pu_var.size(), vector<int>());
    vector<int> temp_vec(n_reads, -1);
    
    for (int i = 0; i < (int)pu_read.size(); ++i){
        // fill in temp_vec
        for (int j = 0; j < (int)pu_read[i].size(); ++j)
            temp_vec[pu_read[i][j]] = i;
        
        // only keep variants covered by reads
        bool is_done = false;
        for (int k = 0; k <= 3; ++k){
            if (4*i+k > pu_var.size()-1){
                is_done = true;
                break;
            }
            for (int j = 0; j < (int)pu_var[4*i+k].size(); ++j){
                if (temp_vec[pu_var[4*i+k][j]] == i)
                    pu_var_ft[4*i+k].push_back(pu_var[4*i+k][j]);
            }
        }
        if (is_done)
            break;
    }
    
    return pu_var_ft;
}

// get consensus
inline void get_consensus(ConsensusSeq &cons, const vector<int> &pu_var_count, const vector<int> &pu_read_count, int start, int end, int min_cvg = 20)
{
    if (floor(double(pu_var_count.size()-1) / 4) > (int)pu_read_count.size() - 1){
        cout << "pu_var_count.size() " << pu_var_count.size() << endl;
        cout << "pu_read_count.size()" << pu_read_count.size() << endl;
        throw runtime_error("ann_clust: floor(double(pu_var_count.size()-1) / 4) > pu_read_count.size() - 1");
    }
    //cons.pu_var_count = pu_var_count;
    //cons.pu_read_count = pu_read_count;
    cons.start = start;
    cons.end = end;
    vector<double> prop(pu_var_count.size(),-1);
    //cons.prop.resize(pu_var_count.size(),-1);
    
    if (pu_var_count.size()==0)
        return;
    
    int start_code = 4*start;
    int end_code = 4*end+3;
    end_code = end_code <= pu_var_count.size()-1 ? end_code : (int)pu_var_count.size()-1;
  
    if (start_code >= pu_var_count.size() || end_code >= pu_var_count.size()){
        cout << "start_code = " << start_code << endl;
        cout << "end_code = " << end_code << endl;
        cout << "pu_var_count.size() = " << pu_var_count.size() << endl;
        throw runtime_error("start_code >= pu_var_count.size() || end_code >= pu_var_count.size()");
    }
    
    if (start_code < 0 || end_code < 0){
        cout << "start_code = " << start_code << endl;
        cout << "end_code = " << end_code << endl;
        throw runtime_error("start_code < 0 || end_code < 0");
    }
    
    for (auto i = start_code; i < end_code; ++i){
        int i_r = int(i/4);
        if (i_r > start && pu_read_count[i_r] > pu_read_count[i_r-1]){
            cout << "pu_read_count = " << pu_read_count << endl;
            cout << "start_code = " << start_code << endl;
            throw runtime_error("coverage should never decrease");
        }
        
        if (pu_read_count[i_r] >= min_cvg){
            //cons.prop[i] = (double) pu_var_count[i] / pu_read_count[i_r];
            prop[i] = (double) pu_var_count[i] / pu_read_count[i_r];
        }else{
            //cons.prop[i] = -1;
            prop[i] = -1;
            cons.end = i_r;
            break;
        }
        //if (cons.prop[i] > 0.5)
        if (prop[i] > 0.5)
            cons.cons_seq.push_back(i);
        
        // to be removed
        /*if (cons.prop[i] > 0.2 && cons.prop[i] < 0.7){
            cout << i << ":" << cons.prop[i] << endl;
            throw runtime_error("cons.prop[i] > 0.2 && cons.prop[i] < 0.7");
        }*/
    }
}

inline void getdepth_from_paf(string paf_file, string cvg_file)
{
    unordered_map<string, vector<int32_t> > cvg_data;
    unordered_map<string, int64_t> contig_len;
    ifstream fs_paf_file; open_infile(fs_paf_file, paf_file);
    int64_t nlines = 0;
    while(true){
        ++nlines;
        if (nlines % 10000 == 0) cout << nlines << endl;
        
        // read line
        string buf;
        getline(fs_paf_file, buf);
        if (fs_paf_file.eof())
            break;
        vector<string> buf_split = split(buf, '\t');
        if (buf_split.size() < 12)
            throw runtime_error("getdepth_from_paf(): buf_split.size() < 12");
        
        // parse
        string chr = buf_split[5];
        int64_t cur_contig_len = stoll(buf_split[6]);
        int64_t tstart = stoll(buf_split[7]);
        int64_t tend = stoll(buf_split[8]);
        
        contig_len[chr] = cur_contig_len;
        
        // calculate depth
        auto it = cvg_data.find(chr);
        if (it == cvg_data.end())
            cvg_data[chr] = vector<int32_t>(cur_contig_len, 0);
        
        it = cvg_data.find(chr);
        for (auto i = tstart; i < tend; ++i){
            ++it->second[i];
        }
        
    }
    cout << nlines << endl;
    fs_paf_file.close();
    
    // print depth
    ofstream fs_cvg_file; open_outfile(fs_cvg_file, cvg_file);
    for (auto & cur_cvg : cvg_data){
        for (auto i = 0; i < cur_cvg.second.size(); ++i){
            if (cur_cvg.second[i] > 0)
                fs_cvg_file << cur_cvg.first << '\t' << contig_len[cur_cvg.first] << '\t' << i+1 << '\t' << cur_cvg.second[i] << endl;
        }
    }
    fs_cvg_file.close();
    
}

#endif /* pileup_h */
