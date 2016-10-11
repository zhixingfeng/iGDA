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
    CondFreq():x(-1),n(-1),p(-1),log_lr(0), log_bf(0), idx(-1){}
    CondFreq(double a_x, double a_n):x(a_x),n(a_n),p(-1),log_lr(0),log_bf(0), idx(-1){}
    CondFreq(double a_x, double a_n, double a_p, double a_log_lr, double a_log_bf, int a_idx):
            x(a_x),n(a_n),p(a_p),log_lr(a_log_lr),log_bf(a_log_bf),idx(a_idx){}
    double x; // number of 1s given the other locus
    double n; // number of 1s at the other locus
    double p; // x/n
    double log_lr; // log likihood ratio
    double log_bf; // log bayes factor
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
    
    // filter out rows with all -1;
    vector<vector<int> > submat_ft;
    for (int i=0; i<(int)submat.size(); i++){
        for (int j=0; j<(int) submat[i].size(); j++){
            if (submat[i][j]!=-1){
                submat_ft.push_back(submat[i]);
                break;
            }
        }
    }
    
    return submat_ft;
}

inline vector<CondFreq> getcondfreq(vector<int> &candidates, vector<vector<int> > &encode_data,
                             vector<ReadRange> &reads_range, vector<int> &pu,
                             vector<double> &p0, string mode="snv")
{
    if (candidates.back() > pu.size())
        throw runtime_error("candidates.back() > pu.size()");
    if (candidates.back() > p0.size())
        throw runtime_error("candidates.back() > p0.size()");
    
    // get marginal frequency for each candidate
    vector<int> freq_mar(candidates.size(), -1);
    for (int i=0; i<(int)candidates.size(); i++)
        freq_mar[i] = pu[candidates[i]-1];
    
    // get p0 for candidates
    vector<double> cand_p0(candidates.size(), -1);
    for (int i=0; i<(int)candidates.size(); i++)
        cand_p0[i] = p0[candidates[i]-1];
    
    // get binary in the subspace spaned by candidates
    vector<vector<int> > submat= getsubspace(candidates, encode_data, (int)pu.size(), reads_range);
    
    // get conditional freq
    vector<CondFreq> condfreq(candidates.size(),CondFreq());
    for (int j=0; j<(int)submat[0].size(); j++){
        // get index of 1 and -1(missing value) at locus j
        vector<int> idx_1;
        vector<int> idx_1n;
        for (int i=0; i<(int)submat.size(); i++){
            if (submat[i][j]==1)
                idx_1.push_back(i);
            if (submat[i][j]==-1)
                idx_1n.push_back(i);
        }
        if (idx_1.size()==0) continue;
        
        // get conditional frquency with locus k
        for (int k=0; k<(int)submat[0].size(); k++){
            if (k==j) continue;
            CondFreq cur_condfreq(0,freq_mar[k]);
            for (int i=0; i<(int) idx_1.size(); i++)
                if (submat[idx_1[i]][k]==1)
                    cur_condfreq.x++;
            for (int i=0; i<(int) idx_1n.size(); i++)
                if (submat[idx_1n[i]][k]==1)
                    cur_condfreq.n--;
            if (cur_condfreq.x > cur_condfreq.n) throw runtime_error("x>n");
            
            cur_condfreq.log_bf = binom_log_bf(cur_condfreq.x, cur_condfreq.n, p0[j]);
            if (cur_condfreq.log_bf > condfreq[j].log_bf){
                condfreq[j] = cur_condfreq;
                condfreq[j].idx = k;
            }
        }
    }
    
    return condfreq;
}


#endif /* condfreq_h */




