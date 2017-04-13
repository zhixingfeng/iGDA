//
//  errormodelsnv.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/10.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef errormodelsnv_h
#define errormodelsnv_h

#include "errormodel.h"

struct BaseFreq{
    BaseFreq(): context(pair<string, string> ("","")), cvg(0), ref('$')
    {
        nvar[0]=0; nvar[1]=0; nvar[2]=0; nvar[3]=0;
    }
    pair<string, string> context; // upstream and downstream context
    int64_t nvar[4]; // 0:A, 1:C, 2:G, 3:T 
    int64_t cvg; // coverage 
    char ref;
};

typedef unordered_map<string, unordered_map<string, vector<BaseFreq> > > ContextEffectAll;
typedef unordered_map<string, unordered_map<string, BaseFreq > > ContextEffect;

class ErrorModelSNV : public ErrorModel {
public:
    
    ErrorModelSNV() : ErrorModel() {}
    virtual ~ErrorModelSNV() {}
    
    void learn(string align_file, string out_prefix);
    void merge(vector<string> &context_files);
        
protected:
    int get_genomesize(string align_file);
    
    void pileup_reads(string align_file, vector<BaseFreq> &pileup);
    
    void print_pileup(string out_file, const vector<BaseFreq> &pileup);

    void get_context_effect_all (const vector<BaseFreq> &pileup, ContextEffectAll &context_effect_all);
    void get_context_effect(const vector<BaseFreq> &pileup, ContextEffect &context_effect);
    
    void print_context_effect_all(string out_file, ContextEffectAll &context_effect_all);
    void print_context_effect(string out_file, ContextEffect &context_effect);

    bool get_context_m5(int locus, int left, int right, const string &tAlignedSeq, pair<string,string> &context);
    
    
protected:
    
};

#endif /* errormodelsnv_hpp */
