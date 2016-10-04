//
//  freqsetminer.cpp
//  iGDA
//
//  Created by Zhixing Feng on 9/26/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "freqsetminer.h"

//void FreqSetMiner::readLineEncode(ifstream & fs_encodefile, vector<int> & encode_line)
//{
    
//}
void FreqSetMiner::readLineCmpReads(ifstream & fs_cmpreads_file, CmpReads & cmpreads_line)
{
    string buf;
    getline(fs_cmpreads_file, buf);
    vector<string> buf_vec = split(buf, '\t');
    if (buf_vec.size() != 4)
        throw runtime_error("error in readLineCmpReads: buf_vec size if not 4.");
    
    cmpreads_line.readsID = split_int(buf_vec[0], ',');
    if (cmpreads_line.readsID.size() != 2)
        throw runtime_error("error in readLineCmpReads: cmpreads_line.readsID size is not 2");
    
    cmpreads_line.cons_code = split_int(buf_vec[1], ',');
    
    cmpreads_line.range_read1 = split_int(buf_vec[2], ',');
    if (cmpreads_line.range_read1.size() != 2)
        throw runtime_error("error in readLineCmpReads: cmpreads_line.range_read1 size is not 2");
    
    cmpreads_line.range_read2 = split_int(buf_vec[3], ',');
    if (cmpreads_line.range_read2.size() != 2)
        throw runtime_error("error in readLineCmpReads: cmpreads_line.range_read2 size is not 2");

    
}

