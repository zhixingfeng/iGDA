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

class DetectSingle
{
    
public:
    DetectSingle(){}
    virtual ~DetectSingle(){}
    
public:
    virtual void loadcontexteffect(string contexteffect_file, int min_context_cvg) = 0;
    virtual void detect(string pileupf_ile, string out_file, double min_prop, int min_cvg) = 0;
    
protected:
    unordered_map<string, unordered_map<string, double> > contexteffect;
};



#endif /* defined(__iGDA__detectsingle__) */
