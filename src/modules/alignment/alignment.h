//
//  alignment.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/11/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef alignment_h
#define alignment_h

#include "../../../include/headers.h"
#include "../../misc/io.h"
#include "../../misc/basic.h"
#include "../tools/tools.h"



class Alignment
{
public:
    Alignment(){}
    virtual ~Alignment(){}
    
public:
    // local_align is a wrapper of Complete-Striped-Smith-Waterman-Library with blasr's (Marc Chaisson's version) default penalties
    inline void local_align(const string & read, const string &ref, StripedSmithWaterman::Alignment &result,
                     int match=5, int mismatch=6, int gap_open=5, int gap_extension=5)
    {
        StripedSmithWaterman::Aligner aligner(match, mismatch, gap_open, gap_extension);
        StripedSmithWaterman::Filter filter;
        StripedSmithWaterman::Alignment alignment;
        
        aligner.Align(read.c_str(), ref.c_str(), (int)ref.size(), filter, &result, 1);

    }

};
#endif /* alignment_h */
