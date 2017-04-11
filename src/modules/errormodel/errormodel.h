//
//  errormodel.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/10.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef errormodel_h
#define errormodel_h

#include "../../../include/headers.h"
#include "../../misc/misc.h"
#include "../alignreader/alignreaderm5.h"
#include "../aligncoder/aligncodersnv.h"

class ErrorModel {
    
public:
    ErrorModel(){}
    virtual ~ErrorModel(){}
    
    virtual void learn(string align_file, string out_prefix)=0;

protected:
    

    
};


#endif /* errormodel_hpp */
