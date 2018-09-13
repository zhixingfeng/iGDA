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
    ErrorModel():left_len(1),right_len(1){}
    virtual ~ErrorModel(){}
    
    virtual inline void set_context_size(int left_len, int right_len)=0;
    virtual void learn(string align_file, string out_prefix)=0;
    virtual void merge(vector<string> &context_files)=0;
    virtual void merge_all(vector<string> &context_all_files)=0;
    
protected:
    
    int left_len;
    int right_len;
    
};


#endif /* errormodel_hpp */
