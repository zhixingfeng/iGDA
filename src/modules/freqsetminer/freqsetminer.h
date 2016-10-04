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
#include "../aligncoder/aligncoder.h"


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
    
    FreqSetMiner(){ptr_aligncoder=NULL;}
    virtual ~FreqSetMiner(){}
    
    void setAlignCoder(AlignCoder *a_ptr_aligncoder){ptr_aligncoder = a_ptr_aligncoder;}
    
    
    virtual vector<int> detectVariantsCoarse(string encode_file, string align_file, string cmpreads_file, double p_cutoff)=0;
    virtual vector<double> detectVariantsSingle(string encode_file, string align_file)=0;
    
    
protected:
    
    //void readLineEncode(ifstream & fs_encodefile, vector<int> & encode_line);
    void readLineCmpReads(ifstream & fs_cmpreads_file, CmpReads & cmpreads_line);

protected:
    AlignCoder *ptr_aligncoder;
};



#endif /* defined(__iGDA__freqsetminer__) */
