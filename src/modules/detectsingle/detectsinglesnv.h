//
//  detectsinglesnv.h
//  iGDA
//
//  Created by Zhixing Feng on 7/5/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__detectsinglesnv__
#define __iGDA__detectsinglesnv__

#include "detectsingle.h"

class DetectSingleSNV : public DetectSingle
{
    
public:
    DetectSingleSNV() : DetectSingle(){}
    virtual ~DetectSingleSNV(){}
    
public:
    void loadcontexteffect(string contexteffect_file, int min_context_cvg = 500);
    void savecontexteffect(string outfile);
    void detect(string pileup_file, string out_file, double min_log_bf = log(50), double min_prop= 0.01, int min_cvg = 20);
    
protected:
    void print_result(string out_file, double min_log_bf = log(50), double min_prop = 0.01, int min_cvg = 20);
};

#endif /* defined(__iGDA__detectsinglesnv__) */
