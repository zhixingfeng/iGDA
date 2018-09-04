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



bool AlignCoderSNV::recode(string m5_file, string var_file, string recode_file, int left_len, int right_len)
{
    // load var_file
    vector<VarData> var_data;
    ifstream fs_varfile;
    int64_t max_code = -1;
    open_infile(fs_varfile, var_file);
    while(true){
        string buf;
        getline(fs_varfile, buf);
        if(fs_varfile.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in " + var_file);
        var_data.push_back(VarData(stod(buf_vec[0]), buf_vec[1][0], stod(buf_vec[2])));
        
        if (stod(buf_vec[2]) > max_code)
            max_code = stod(buf_vec[2]);
    }
    fs_varfile.close();
    
    // fill template of var_data
    vector<bool> var_data_temp(max_code + 1, false);
    for (int64_t i = 0; i < var_data.size(); ++i)
        var_data_temp[var_data[i].code] = true;
    
    
    // scan m5_file and recode
    if (p_alignreader==NULL)
        throw runtime_error("AlignCoderSNV::recode(): p_alignreader has not be set.");
    ofstream p_outfile;
    open_outfile(p_outfile, recode_file);
    
    p_alignreader->open(m5_file);
    Align align;
    int nline = 0;
    while(p_alignreader->readline(align)){
        ++nline;
        //cout << nline << endl;
        
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
            //cout << "i=" << i <<',';
            if (align.tAlignedSeq[i]=='-')
                continue;
            
            if (4*cur_pos+3 > max_code)
                break;
            // realign if hit detected variants
            int score_A = -1000000;
            int score_C = -1000000;
            int score_G = -1000000;
            int score_T = -1000000;
            int score_ref = -1000000;
            bool is_var = false;
            
            string cur_qseq;
            string cur_rseq;
            pair<string, string> context;
            if (var_data_temp[4*cur_pos] || var_data_temp[4*cur_pos+1] || var_data_temp[4*cur_pos+2] || var_data_temp[4*cur_pos+3]){
                bool rl = this->get_context_m5(i, left_len, right_len, align.tAlignedSeq, context);
                if (!rl){
                    ++cur_pos;
                    continue;
                }
                
                for (auto j = i - context.first.size() + 1; j <= i + context.second.size(); ++j){
                    if (align.qAlignedSeq[j]!='-')
                        cur_qseq.push_back(align.qAlignedSeq[j]);
                }
                string cur_rseq = context.first + context.second;
                seqan::Align<string, seqan::ArrayGaps> cur_realign;
                score_ref = this->realign(cur_realign, cur_qseq, cur_rseq);
                
            }
            
            if (var_data_temp[4*cur_pos]){
                if (align.tAlignedSeq[i] == 'A'){
                    char err_msg[1000];
                    sprintf(err_msg, "read %d: detect variant is equal to reference. var is A, and ref at %dth alignment is %c at locus %d", nline, i, align.tAlignedSeq[i], cur_pos);
                    throw runtime_error(err_msg);
                }
                is_var = true;
            }
            
            if (var_data_temp[4*cur_pos+1]){
                if (align.tAlignedSeq[i] == 'C'){
                    char err_msg[1000];
                    sprintf(err_msg, "read %d: detect variant is equal to reference. var is C, and ref at %dth alignment is %c at locus %d", nline, i, align.tAlignedSeq[i], cur_pos);
                    throw runtime_error(err_msg);
                }
                is_var = true;
            }
            
            if (var_data_temp[4*cur_pos+2]){
                if (align.tAlignedSeq[i] == 'G'){
                    char err_msg[1000];
                    sprintf(err_msg, "read %d: detect variant is equal to reference. var is G, and ref at %dth alignment is %c at locus %d", nline, i, align.tAlignedSeq[i], cur_pos);
                    throw runtime_error(err_msg);
                }
                is_var = true;
            }
            
            if (var_data_temp[4*cur_pos+3]){
                if (align.tAlignedSeq[i] == 'T'){
                    char err_msg[1000];
                    sprintf(err_msg, "read %d: detect variant is equal to reference. var is T, and ref at %dth alignment is %c at locus %d", nline, i, align.tAlignedSeq[i], cur_pos);
                    throw runtime_error(err_msg);
                }
                is_var = true;
            }
            
            ++cur_pos;
        }
        p_outfile << endl;
    }
    
    p_alignreader->close();
    
    
    p_outfile.close();
    
    
    return true;
}

bool AlignCoderSNV::get_context_m5(int i, int left, int right, const string &tAlignedSeq, pair<string,string> &context)
{
    int n_base, k;
    
    // search left
    n_base = 0; k = 0;
    char cur_base = tAlignedSeq[i];
    string context_left("");
    while (n_base<=left){
        k++;
        if (i - k < 0)
            return false;
        if (tAlignedSeq[i-k] != cur_base && tAlignedSeq[i-k] != '-'){
            cur_base = tAlignedSeq[i-k];
            n_base++;
        }
        
    }
    k--;
    while (k>=0){
        if (tAlignedSeq[i-k]!='-')
            context_left.push_back(tAlignedSeq[i-k]);
        k--;
    }
    
    // search right
    n_base = 0; k = 0;
    cur_base = tAlignedSeq[i];
    string context_right("");
    while (n_base <= right){
        k++;
        if (i + k >= tAlignedSeq.size())
            return false;
        if (tAlignedSeq[i+k] != cur_base && tAlignedSeq[i+k] != '-'){
            cur_base = tAlignedSeq[i+k];
            n_base++;
        }
    }
    k--;
    for (int j=i+1; j<=i+k; j++){
        if (tAlignedSeq[j]!='-')
            context_right.push_back(tAlignedSeq[j]);
    }
    context.first = context_left;
    context.second = context_right;
    return true;
}


int AlignCoderSNV::realign(seqan::Align<string, seqan::ArrayGaps> &cur_realign, const string &qseq, const string &rseq)
{
    seqan::resize(seqan::rows(cur_realign),2);
    seqan::assignSource(seqan::row(cur_realign, 0), qseq);
    seqan::assignSource(seqan::row(cur_realign, 1), rseq);
    
    int cur_score = seqan::globalAlignment(cur_realign, seqan::Score<int, seqan::Simple>(2, -4, -2, -4), seqan::AffineGaps());
    return cur_score;
}








