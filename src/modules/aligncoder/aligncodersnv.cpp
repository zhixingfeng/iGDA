//
//  aligncodersnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "aligncodersnv.h"
#include "../alignreader/alignreaderm5.h"

bool AlignCoderSNV::encode(string alignfile, string outfile)
{
    ofstream p_outfile;
    open_outfile(p_outfile, outfile);
    AlignReaderM5 alignreaderm5;
    AlignReader *p_alignreader = &alignreaderm5;
    
    p_alignreader->open(alignfile);
    Align align;
    int nline = 0;
    while(p_alignreader->readline(align)){
        ++nline;
        int alen = (int) align.matchPattern.size();
        if ( !(align.qAlignedSeq.size()==alen && align.tAlignedSeq.size()==alen) )
            throw runtime_error("incorrect match patter in line " + to_string(nline));
        
        int cur_pos = align.tStart;
        for (int i=0; i<alen; i++){
            if (align.tAlignedSeq[i]!=align.qAlignedSeq[i] && align.tAlignedSeq[i]!='-' && align.qAlignedSeq[i]!='-')
                p_outfile << cur_pos << '\t';
            if (align.tAlignedSeq[i]!='-')
                cur_pos++;
        }
        p_outfile << endl;
    }
    
    p_alignreader->close();
    
    
    p_outfile.close();
    return true;
}


