//
//  aligncodersnv.h
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__aligncodersnv__
#define __iGDA__aligncodersnv__

#include "aligncoder.h"

class AlignCoderSNV : public AlignCoder
{
public:
    AlignCoderSNV() : AlignCoder() {};
    virtual ~AlignCoderSNV(){};
    
    bool encode(string alignfile, string outfile);
    pair<int, char> decode(int code);
    
    // code formula: shift is A=+0, C=+1, G=+2, T=+3, code=4*(pos-1) + shift
    int binary_code(int pos, char base);
    
    // encode alignment result of ssw, start_pos is 0-based !
    bool encode(const StripedSmithWaterman::Alignment &alignment, const string &read, const string &ref,
                    int start_pos, vector<int> &encode_data);

    // recode according to detected variants
    bool recode(string m5_file, string var_file, string recode_file, int left_len, int right_len);
    
protected:
    bool get_context_m5(int locus, int left, int right, const string &tAlignedSeq, pair<string,string> &context);
    int realign(seqan::Align<string, seqan::ArrayGaps> &cur_realign, const string &qseq, const string &rseq);
};

#endif /* defined(__iGDA__aligncodersnv__) */
