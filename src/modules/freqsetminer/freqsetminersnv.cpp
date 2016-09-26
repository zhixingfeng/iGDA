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
    ifstream fs_contigfile; open_infile(fs_contigfile, cmpreadsfile);
    
    
    
    fs_encodefile.close();
    fs_contigfile.close();
    
    return true;

}