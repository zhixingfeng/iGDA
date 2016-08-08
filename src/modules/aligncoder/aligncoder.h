//
//  aligncoder.h
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__aligncoder__
#define __iGDA__aligncoder__

#include <headers.h>

class AlignCoder
{
public:
    AlignCoder(){};
    virtual ~AlignCoder(){};
    
    virtual bool encode(string alignfile, string outfile)=0;
    virtual pair<int, char> decode(int code)=0; // int is position, char is base
};


#endif /* defined(__iGDA__aligncoder__) */
