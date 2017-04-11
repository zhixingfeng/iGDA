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
    BaseFreq(): context(pair<string, string> ("","")), cvg(0)
    {
        nvar[0]=0; nvar[1]=0; nvar[2]=0; nvar[3]=0;
    }
    pair<string, string> context; // upstream and downstream context
    int nvar[4]; // 0:A, 1:C, 2:G, 3:T 
    int cvg; // coverage 
};


class ErrorModelSNV : public ErrorModel {
public:
    
    ErrorModelSNV() : ErrorModel() {}
    virtual ~ErrorModelSNV() {}
    
    void learn(string align_file, string out_prefix);
        
protected:
    int get_genomesize(string align_file);
    
    void pileup_reads(string align_file, vector<BaseFreq> &pileup);
    
    void print_pileup(string out_file, const vector<BaseFreq> &pileup);

    
    bool get_context_m5(int locus, int left, int right, const string &tAlignedSeq, pair<string,string> &context);
    
    
protected:
    
};

#endif /* errormodelsnv_hpp */
