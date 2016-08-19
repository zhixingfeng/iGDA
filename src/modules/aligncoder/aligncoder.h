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
#include "../alignreader/alignreader.h"
#include "../alignreader/alignreaderm5.h"

class AlignCoder
{
public:
    AlignCoder(){p_alignreader = NULL;}
    virtual ~AlignCoder(){};
    
    virtual bool encode(string alignfile, string outfile)=0;
    virtual pair<int, char> decode(int code)=0; // int is position, char is base
    
    inline void setAlignReader(AlignReader * a_p_alignreader){p_alignreader = a_p_alignreader;}
    
protected:
    AlignReader *p_alignreader;
};


#endif /* defined(__iGDA__aligncoder__) */
