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
    AlignCoderSNV() : AlignCoder() {};
    virtual ~AlignCoderSNV(){};
    
    bool encode(string alignfile, string outfile);
    pair<int, char> decode(int code);
    
    // code formula: shift is A=+0, C=+1, G=+2, T=+3, code=4*(pos-1) + shift
    int binary_code(int pos, char base); 

};

#endif /* defined(__iGDA__aligncodersnv__) */
