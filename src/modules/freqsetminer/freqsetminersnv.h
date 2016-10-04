//
//  freqsetminersnv.h
//  iGDA
//
//  Created by Zhixing Feng on 9/26/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__freqsetminersnv__
#define __iGDA__freqsetminersnv__

#include "../../misc/misc.h"
#include "freqsetminer.h"

class FreqSetMinerSNV : public FreqSetMiner
{
public:
    FreqSetMinerSNV() : FreqSetMiner() {}
    virtual ~FreqSetMinerSNV(){}
    
    vector<int> detectVariantsCoarse(string encode_file, string align_file, string cmpreads_file, double p_cutoff);
    
    vector<double> detectVariantsSingle(string encode_file, string align_file);
};

#endif /* defined(__iGDA__freqsetminersnv__) */
