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
                p_outfile << this->binary_code(cur_pos, align.qAlignedSeq[i]) << '\t';
            if (align.tAlignedSeq[i]!='-')
                cur_pos++;
        }
        p_outfile << endl;
    }
    
    p_alignreader->close();
    
    
    p_outfile.close();
    return true;
}

pair<int, char> AlignCoderSNV::decode(int code)
{
    pair<int, char> rl;
    rl.first = int ((code-1) / 4) + 1;
    int shift = code % 4;
    switch(shift){
        case 1:
            rl.second = 'A';
            break;
        case 2:
            rl.second = 'C';
            break;
        case 3:
            rl.second = 'G';
            break;
        case 0:
            rl.second = 'T';
            break;
        default:
            throw runtime_error("AlignCoderSNV::decode: incorrect code");
    }
    return rl;
}
int AlignCoderSNV::binary_code(int pos, char base)
{
    int shift = 0;
    switch(base){
        case 'A':
            shift = 1;
            break;
        case 'C':
            shift = 2;
            break;
        case 'G':
            shift = 3;
            break;
        case 'T':
            shift = 4;
            break;
        default:
            throw runtime_error("AlignCoderSNV::binary_code: incorrect base");
    }
    
    return 4*(pos-1) + shift;
}




