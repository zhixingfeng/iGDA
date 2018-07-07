//
//  detectsingle.h
//  iGDA
//
//  Created by Zhixing Feng on 7/5/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__detectsingle__
#define __iGDA__detectsingle__
#include "../errormodel/errormodel.h"
#include "../errormodel/errormodelsnv.h"

struct DetectSingleResult
{
    // orignal data of pileup
    int64_t locus;
    string context_left;
    string context_right;
    char ref_base;
    int A_count, C_count, G_count, T_count;
    int cvg;
    
    // detection
    double A_freq, C_freq, G_freq, T_freq;
    double A_bf, C_bf, G_bf, T_bf;
    double A_freq_ref, C_freq_ref, G_freq_ref, T_freq_ref, cvg_ref;
    
    
};

class DetectSingle
{
    
public:
    DetectSingle(){}
    virtual ~DetectSingle(){}
    
public:
    virtual void loadcontexteffect(string contexteffect_file, int min_context_cvg = 500) = 0;
    virtual void savecontexteffect(string outfile) = 0;
    virtual void detect(string pileup_file, string out_file, double min_bf = 10, double min_prop= 0.02, int min_cvg = 20) = 0;
    
protected:
    unordered_map<string, unordered_map<string, vector<double> > > contexteffect;
    vector<DetectSingleResult> result;
};



#endif /* defined(__iGDA__detectsingle__) */
