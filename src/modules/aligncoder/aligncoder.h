//
//  aligncoder.h
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__aligncoder__
#define __iGDA__aligncoder__

#include <headers.h>
#include "../alignreader/alignreader.h"
#include "../alignreader/alignreaderm5.h"
#include "../../../tools/tools.h"
#include <seqan/align.h>
#include "../../misc/io.h"

struct VarData
{
public:
    VarData(int64_t a_locus, int64_t a_base, int64_t a_code):locus(a_locus), base(a_base), code(a_code){}
    int64_t locus;
    char base;
    int64_t code;
};


class AlignCoder
{
public:
    AlignCoder(){p_alignreader = NULL;}
    virtual ~AlignCoder(){};
    
    virtual bool encode(string alignfile, string outfile)=0;
    virtual pair<int, char> decode(int code)=0; // int is position, char is base
    virtual int binary_code(int pos, char base)=0; 
    
    virtual bool recode_legacy(string m5_file, string var_file, string recode_file, int left_len, int right_len, bool is_report_ref = true) = 0;
    virtual bool recode(string m5_file, string var_file, string recode_file, int left_len, int right_len, bool is_report_ref = true) = 0;
    inline void setAlignReader(AlignReader * a_p_alignreader){p_alignreader = a_p_alignreader;}
    
protected:
    AlignReader *p_alignreader;
};


#endif /* defined(__iGDA__aligncoder__) */
