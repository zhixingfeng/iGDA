//
//  aligncodersnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "aligncodersnv.h"
#include "../alignreader/alignreader.h"
#include "../alignreader/alignreaderm5.h"
#include "../../misc/misc.h"

// public

bool AlignCoderSNV::encode(string alignfile, string outfile)
{
    if (p_alignreader==NULL)
        throw runtime_error("AlignCoderSNV::encode: p_alignreader has not be set.");
    ofstream p_outfile;
    open_outfile(p_outfile, outfile);
    
    p_alignreader->open(alignfile);
    Align align;
    int nline = 0;
    while(p_alignreader->readline(align)){
        ++nline;
        
        // expections
        int alen = (int) align.matchPattern.size();
        if ( !(align.qAlignedSeq.size()==alen && align.tAlignedSeq.size()==alen) )
            throw runtime_error("incorrect match patter in line " + to_string(nline));
        if (align.qStrand != '+')
            throw runtime_error("qStrand should be + in line " + to_string(nline));
        
        // reverse alignment if it is aligned to negative strand
        if (align.tStrand != '+'){
            align.qAlignedSeq = getrevcomp(align.qAlignedSeq);
            align.tAlignedSeq = getrevcomp(align.tAlignedSeq);
        }
        
        // encode
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
    rl.first = int (code / 4);
    int shift = code % 4;
    switch(shift){
        case 0:
            rl.second = 'A';
            break;
        case 1:
            rl.second = 'C';
            break;
        case 2:
            rl.second = 'G';
            break;
        case 3:
            rl.second = 'T';
            break;
        default:
            throw runtime_error("AlignCoderSNV::decode: incorrect code");
    }
    return rl;
}

// private
int AlignCoderSNV::binary_code(int pos, char base)
{
    int shift = 0;
    switch(base){
        case 'A':
            shift = 0;
            break;
        case 'C':
            shift = 1;
            break;
        case 'G':
            shift = 2;
            break;
        case 'T':
            shift = 3;
            break;
        default:
            throw runtime_error("AlignCoderSNV::binary_code: incorrect base");
    }
    
    return 4*pos + shift;
}

bool AlignCoderSNV::encode(const StripedSmithWaterman::Alignment &alignment, const string &read, const string &ref,
                               int start_pos, vector<int> &encode_data)
{
    // scan cigar to encode the variants (only work for positive strand !!)
    int read_pos = alignment.query_begin;
    int ref_pos = alignment.ref_begin;
    for (int i=0; i<(int)alignment.cigar.size(); ++i){
        // parse cigar (high 28 bits: length, low 4 bits: M/I/D/S/X (0/1/2/4/8)), = is 7
        uint32_t cigar_type = alignment.cigar[i] & 0xF;
        uint32_t cigar_len = alignment.cigar[i] >> 4;
        
        switch(cigar_type){
            // if cigar_type is M, shift read_pos and ref_pos
            case 0:
                read_pos += cigar_len;
                ref_pos += cigar_len;
                break;
            // if cigar_type is I, shift read_pos
            case 1:
                read_pos += cigar_len;
                break;
            // if cigar_type is D, shift ref_pos
            case 2:
                ref_pos += cigar_len;
                break;
            // if cigar_type is S, do nothing
            case 4:
                break;
            // if cigar_type is =, shift read_pos and ref_pos    
            case 7:
                read_pos += cigar_len;
                ref_pos += cigar_len;
                break;

            // if cigar_type is X, encode the variant and shift read_pos and ref_pos
            case 8:
                for (int j=0; j<(int)cigar_len; ++j){
                    encode_data.push_back(binary_code(start_pos + ref_pos + j, read[read_pos+j]));
                }
                read_pos += cigar_len;
                ref_pos += cigar_len;
                break;
            default:
                throw runtime_error("unknow cigar type.");
                break;
            
        }
    }
    return true;
}


