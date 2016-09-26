//
//  freqsetminer.h
//  iGDA
//
//  Created by Zhixing Feng on 9/26/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__freqsetminer__
#define __iGDA__freqsetminer__

#include <headers.h>
#include "../../../tools/tools.h"


struct CmpReads
{
    pair<int, int> readsID;
    vector<int> cons_code;
};

class FreqSetMiner
{
public:
    FreqSetMiner(){}
    virtual ~FreqSetMiner(){}
    
    virtual bool mapEncodetoCmpReads(string encodefile, string cmpreadsfile)=0;

    void readLineEncode(const ifstream & fs_encodefile, vector<int> & encode_line);
    void readLineCmpReads(const ifstream & fs_cmpreadsfile, CmpReads & cmpreads_line);
};



#endif /* defined(__iGDA__freqsetminer__) */
