//
//  freqsetminersnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 9/26/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "freqsetminersnv.h"

//naive hash table based approach
bool FreqSetMinerSNV::mapEncodetoCmpReads(string encodefile, string cmpreadsfile)
{
    ifstream fs_encodefile; open_infile(fs_encodefile, encodefile);
    ifstream fs_cmpreadsfile; open_infile(fs_cmpreadsfile, cmpreadsfile);
    
    //vector<int> encode_line;
    //readLineEncode(fs_encodefile, encode_line);
    CmpReads cmpreads_line;
    readLineCmpReads(fs_cmpreadsfile, cmpreads_line);
    
    fs_encodefile.close();
    fs_cmpreadsfile.close();
    
    return true;

}
