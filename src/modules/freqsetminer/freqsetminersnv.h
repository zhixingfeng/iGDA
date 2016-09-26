//
//  freqsetminersnv.h
//  iGDA
//
//  Created by Zhixing Feng on 9/26/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__freqsetminersnv__
#define __iGDA__freqsetminersnv__

#include "freqsetminer.h"

class FreqSetMinerSNV : public FreqSetMiner
{
public:
    FreqSetMinerSNV() : FreqSetMiner() {}
    virtual ~FreqSetMinerSNV(){}
    
    bool mapEncodetoCmpReads(string encodefile, string cmpreadsfile);
};

#endif /* defined(__iGDA__freqsetminersnv__) */
