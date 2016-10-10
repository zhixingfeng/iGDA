//
//  condfreq.h
//  iGDA
//
//  Created by Zhixing Feng on 16/10/9.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef condfreq_h
#define condfreq_h
#include "io.h"
#include "../modules/aligncoder/aligncodersnv.h"

struct CondFreq
{
    CondFreq():x(-1),n(-1),p(-1),log_lr(0),idx(-1){}
    CondFreq(double a_x, double a_n):x(a_x),n(a_n),p(-1),log_lr(0),idx(-1){}
    CondFreq(double a_x, double a_n, double a_p, double a_log_lr, int a_idx):x(a_x),n(a_n),p(a_p),log_lr(a_log_lr),idx(a_idx){}
    double x; // number of 1s given the other locus
    double n; // number of 1s at the other locus
    double p; // x/n
    double log_lr; // likihood ratio
    int idx; // id of the other locus that give the most significant result
    
};

inline vector<vector<int> > getsubspace(vector<int> &candidates, vector<vector<int> > &encode_data, int max_encode,
                                        vector<ReadRange> &reads_range, string mode="snv")
{
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data != reads_range.");
    
    // setup encoder
    AlignCoder *aligncoder=NULL;
    AlignCoderSNV aligncodersnv;
    
    if (mode == "snv") aligncoder = &aligncodersnv;
    if (aligncoder == NULL)
        throw runtime_error("illegal mode.");
    
    // get loci of candidates
    vector<int> cand_loci(candidates.size(),-1);
    for (int i=0; i<(int)candidates.size(); i++)
        cand_loci[i] = aligncoder->decode(candidates[i]).first;
    
    // initialize the binary matrix in the subspace
    vector<vector<int> > submat;
    for (int i=0; i<(int)reads_range.size(); i++){
        vector<int> cur_submat(candidates.size(),0);
        for (int j=0; j<(int)candidates.size(); j++)
            if (cand_loci[j] < reads_range[i].first || cand_loci[j] > reads_range[i].second)
                cur_submat[j] = -1;
        submat.push_back(cur_submat);
    }
    
    // create a index vector for candidates
    vector<int> index(max_encode, -1);
    for (int i=0; i<(int)candidates.size(); i++)
        index[candidates[i]-1] = i;
    
    
    // create the binary matrix in the subspace
    for (int i=0; i<(int)encode_data.size(); i++){
        for (int j=0; j<(int)encode_data[i].size(); j++){
            if (index[encode_data[i][j]-1] != -1){
                int k = index[encode_data[i][j]-1];
                submat[i][k] = 1;
            }
        }
    }
    
    return submat;
}

/*inline vector<CondFreq> getcondfreq(vector<int> &candidates, vector<vector<int> > &encode_data,
                                    vector<ReadRange> &reads_range, vector<int> &pu,
                                    vector<double> &p0, string mode="snv")
{
    if (candidates.back() > pu.size())
        throw runtime_error("candidates.back() > pu.size()");
    
    // get marginal frequency for each candidate
    vector<int> freq_mar(candidates.size(), -1);
    for (int i=0; i<(int)candidates.size(); i++)
        freq_mar[i] = pu[candidates[i]-1];
    
    //
    vector<CondFreq> condfreq;
    
    
    
    return condfreq;
}*/

#endif /* condfreq_h */




