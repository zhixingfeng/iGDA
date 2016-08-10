//
//  m5tofa.h
//  iGDA
//
//  Created by Zhixing Feng on 8/10/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef iGDA_m5tofa_h
#define iGDA_m5tofa_h

#include "../modules/modules.h"
#include "../../include/headers.h"

inline bool m5tofa(string m5file, string fafile)
{
    AlignReaderM5 alignreaderm5;
    Align align;
    
    ofstream p_fafile; open_outfile(p_fafile, fafile);
    alignreaderm5.open(m5file);
    
    while(alignreaderm5.readline(align)) {
        p_fafile << '>' << align.qName << endl;
        string qSeq = "";
        for (int i=0; i<(int)align.qAlignedSeq.size(); i++)
            if (align.qAlignedSeq[i] != '-')
                qSeq.push_back(align.qAlignedSeq[i]);
        p_fafile << qSeq << endl;
    }
    alignreaderm5.close();
    p_fafile.close();
    return true;
}

#endif
