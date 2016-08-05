//
//  aligncodersnv.h
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__aligncodersnv__
#define __iGDA__aligncodersnv__

#include "aligncoder.h"

class AlignCoderSNV : public AlignCoder
{
public:
    AlignCoderSNV(){};
    virtual ~AlignCoderSNV(){};
    
    bool encode(string alignfile, string outfile);

};

#endif /* defined(__iGDA__aligncodersnv__) */
