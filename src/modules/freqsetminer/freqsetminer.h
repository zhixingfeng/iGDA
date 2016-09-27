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
    vector<int> readsID;
    vector<int> cons_code;
    vector<int> range_read1;
    vector<int> range_read2;
};

class FreqSetMiner
{
    
public:
    
    FreqSetMiner(){}
    virtual ~FreqSetMiner(){}
    
    virtual bool mapEncodetoCmpReads(string encodefile, string cmpreadsfile)=0;

protected:
    
    void readLineEncode(ifstream & fs_encodefile, vector<int> & encode_line);
    void readLineCmpReads(ifstream & fs_cmpreadsfile, CmpReads & cmpreads_line);
};



#endif /* defined(__iGDA__freqsetminer__) */
