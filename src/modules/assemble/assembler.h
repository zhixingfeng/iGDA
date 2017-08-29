//
//  assemble.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef assemble_h
#define assemble_h

#include "../../../include/headers.h"
#include "../aligncoder/aligncodersnv.h"

class Assembler 
{
public:
    Assembler(){}
    virtual ~Assembler(){}
    
public:
    void get_variants(string dforest_file, string out_file, double min_condprob);
    
    
};

#endif /* assemble_h */
