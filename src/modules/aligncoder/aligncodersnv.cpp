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
    
    //ofstream p_reffile;
    //open_outfile(p_reffile, outfile + ".ref");
    
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
        /*if (align.tStrand != '+'){
            align.qAlignedSeq = getrevcomp(align.qAlignedSeq);
            align.tAlignedSeq = getrevcomp(align.tAlignedSeq);
        }*/
        
        // encode
        int cur_pos = align.tStart;
        for (int i=0; i<alen; i++){
            if (align.tAlignedSeq[i] == 'A' || align.tAlignedSeq[i] == 'C' || align.tAlignedSeq[i] == 'G' || align.tAlignedSeq[i] == 'T'){
                if (align.qAlignedSeq[i] != '-'){
                    if (align.tAlignedSeq[i] != align.qAlignedSeq[i]){
                        if (align.qAlignedSeq[i] == 'A' || align.qAlignedSeq[i] == 'C' || align.qAlignedSeq[i] == 'G' || align.qAlignedSeq[i] == 'T')
                            p_outfile << this->binary_code(cur_pos, align.qAlignedSeq[i]) << '\t';
                    }else{
                        //p_reffile << this->binary_code(cur_pos, align.tAlignedSeq[i]) << '\t';
                    }
                }
                ++cur_pos;
            }else{
                if (align.tAlignedSeq[i]!='-')
                    ++cur_pos;
            }
            
            /*if (align.tAlignedSeq[i]!=align.qAlignedSeq[i] && align.tAlignedSeq[i]!='-' && align.qAlignedSeq[i]!='-')
                p_outfile << this->binary_code(cur_pos, align.qAlignedSeq[i]) << '\t';
            if (align.tAlignedSeq[i]!='-')
                cur_pos++;*/
        }
        p_outfile << endl;
        //p_reffile << endl;
    }
    
    p_alignreader->close();
    
    
    p_outfile.close();
    //p_reffile.close();
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



bool AlignCoderSNV::recode_legacy(string m5_file, string var_file, string recode_file, int left_len, int right_len, bool is_report_ref)
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
    vector<bool> var_data_temp(max_code + 4, false);
    for (int64_t i = 0; i < var_data.size(); ++i)
        var_data_temp[var_data[i].code] = true;
    
    
    // scan m5_file and recode
    if (p_alignreader==NULL)
        throw runtime_error("AlignCoderSNV::recode(): p_alignreader has not be set.");
    ofstream p_outfile;
    open_outfile(p_outfile, recode_file);
    ofstream p_outfile_ref;
    open_outfile(p_outfile_ref, recode_file + ".ref");
    
    p_alignreader->open(m5_file);
    Align align;
    int nline = 0;
    while(p_alignreader->readline(align)){
        ++nline;
        if (nline % 100 == 0)
            cout << nline << endl;
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
            //cout << "nline=" << nline << ", i=" <<i << endl;
            if (align.tAlignedSeq[i]=='-')
                continue;
            
            if (4*cur_pos+3 > max_code + 3)
                break;
            
            // realign if hit detected variants
            int score_A = MIN_SCORE;
            int score_C = MIN_SCORE;
            int score_G = MIN_SCORE;
            int score_T = MIN_SCORE;
            int score_ref = MIN_SCORE;
            bool is_var = false;
            
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_A;
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_C;
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_G;
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_T;
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_ref;
            
            string cur_qseq;
            string cur_rseq;
            pair<string, string> context;
            char ref_base;
            
            // if the current locus hit a detected variants then realign
            if (var_data_temp[4*cur_pos] || var_data_temp[4*cur_pos+1] || var_data_temp[4*cur_pos+2] || var_data_temp[4*cur_pos+3]){
                bool rl = this->get_context_m5(i, left_len, right_len, align.tAlignedSeq, context);
                if (!rl){
                    ++cur_pos;
                    continue;
                }
                
                is_var = true;
                
                // get left query sequence length
                int64_t k = 0;
                int64_t cur_qseq_start = i;
                while(true){
                    if (align.tAlignedSeq[cur_qseq_start]!='-')
                        k++;
                    if (k >= context.first.size())
                        break;
                    --cur_qseq_start;
                }
                
                
                // get right query sequence length
                k = 0;
                int64_t cur_qseq_end = i+1;
                while(true){
                    if (align.tAlignedSeq[cur_qseq_end]!='-')
                        k++;
                    if (k >= context.second.size())
                        break;
                    ++cur_qseq_end;
                }
                
                if (cur_qseq_start < 0)
                    throw runtime_error("cur_qseq_start < 0");
                if (cur_qseq_end >= alen)
                    throw runtime_error("cur_qseq_end >= alen");
                
                for (auto j = cur_qseq_start; j <= cur_qseq_end; ++j){
                    if (align.qAlignedSeq[j]!='-')
                        cur_qseq.push_back(align.qAlignedSeq[j]);
                }
                
                cur_rseq = context.first + context.second;
                
                if (cur_qseq == ""){
                    ++cur_pos;
                    continue;
                }
                
                if (cur_rseq == "")
                    throw runtime_error("cur_rseq is empty");
                
                ref_base = cur_rseq[context.first.size()-1];
                
                // realign to the reference
                score_ref = this->realign(cur_realign_ref, cur_qseq, cur_rseq);
                
                // realign to variants
                // realign to A
                if (var_data_temp[4*cur_pos]){
                    if (ref_base == 'A')
                        throw runtime_error("reference is A but variants is also A at line " + to_string(nline) + ":" + to_string(i));
                    cur_rseq[context.first.size()-1] = 'A';
                    score_A = this->realign(cur_realign_A, cur_qseq, cur_rseq);
                }
                
                // realign to C
                if (var_data_temp[4*cur_pos + 1]){
                    if (ref_base == 'C')
                        throw runtime_error("reference is C but variants is also C at line " + to_string(nline) + ":" + to_string(i));
                    cur_rseq[context.first.size()-1] = 'C';
                    score_C = this->realign(cur_realign_C, cur_qseq, cur_rseq);
                }
                
                // realign to G
                if (var_data_temp[4*cur_pos + 2]){
                    if (ref_base == 'G')
                        throw runtime_error("reference is G but variants is also G at line " + to_string(nline) + ":" + to_string(i));
                    cur_rseq[context.first.size()-1] = 'G';
                    score_G = this->realign(cur_realign_G, cur_qseq, cur_rseq);
                }
                
                // realign to T
                if (var_data_temp[4*cur_pos + 3]){
                    if (ref_base == 'T')
                        throw runtime_error("reference is T but variants is also T at line " + to_string(nline) + ":" + to_string(i));
                    cur_rseq[context.first.size()-1] = 'T';
                    score_T = this->realign(cur_realign_T, cur_qseq, cur_rseq);
                }
                
            }else{
                ++cur_pos;
                continue;
            }
    
            // to be removed
            /*if (cur_pos == 2827007){
                cout << "read: " << nline << endl;
                cout << "score: " << endl;
                cout << "A = " << score_A << ", " << "C = " << score_C << ", " << "G = " << score_G << ", " << "T = " << score_T << ", ref = " << score_ref << endl;
                
                if (score_A > MIN_SCORE){
                    cout << "realign A" << endl;
                    cout << cur_realign_A << endl;
                }
                
                if (score_C > MIN_SCORE){
                    cout << "realign C" << endl;
                    cout << cur_realign_C << endl;
                }
                if (score_G > MIN_SCORE){
                    cout << "realign G" << endl;
                    cout << cur_realign_G << endl;
                }
                if (score_T > MIN_SCORE){
                    cout << "realign T" << endl;
                    cout << cur_realign_T << endl;
                }
                    
                if (score_ref > MIN_SCORE){
                    cout << "realign ref" << endl;
                    cout << cur_realign_ref << endl;
                }
                cout << "ref_base = " << ref_base << endl;
                getchar();
             }*/
            
            // recode
            if (score_A == MIN_SCORE && score_C == MIN_SCORE && score_G == MIN_SCORE && score_T == MIN_SCORE)
                throw runtime_error("A,C,G,T == MIN_SCORE, no alignment was done");
            
            // A
            if (score_A > score_C && score_A > score_G && score_A > score_T && score_A > score_ref){
                p_outfile << 4*cur_pos << '\t';
            }
            
            // C
            if (score_C > score_A && score_C > score_G && score_C > score_T && score_C > score_ref){
                p_outfile << 4*cur_pos+1 << '\t';
            }
            
            // G
            if (score_G > score_A && score_G > score_C && score_G > score_T && score_G > score_ref){
                p_outfile << 4*cur_pos+2 << '\t';
            }
            
            // T
            if (score_T > score_A && score_T > score_C && score_T > score_G && score_T > score_ref){
                p_outfile << 4*cur_pos+3 << '\t';
            }
            
            // ref
            if (is_report_ref){
                if (score_ref > score_A && score_ref > score_C && score_ref > score_T){
                    switch (ref_base) {
                        case 'A':
                            p_outfile_ref << 4*cur_pos << '\t';
                            break;
                        case 'C':
                            p_outfile_ref << 4*cur_pos + 1 << '\t';
                            break;
                        case 'G':
                            p_outfile_ref << 4*cur_pos + 2 << '\t';
                            break;
                        case 'T':
                            p_outfile_ref << 4*cur_pos + 3 << '\t';
                            break;
                        default:
                            throw runtime_error("ref_base is not A, C, G or T");
                            break;
                    }
                }
            }
            
            ++cur_pos;
        }
        p_outfile << endl;
        p_outfile_ref << endl;
    }
    
    p_alignreader->close();
    
    
    p_outfile.close();
    p_outfile_ref.close();
    
    return true;
}
bool AlignCoderSNV::recode(string m5_file, string var_file, string recode_file, int left_len, int right_len, bool is_report_ref)
{
    AlignReaderM5 alignreaderm5;
    AlignReader *p_alignreader = &alignreaderm5;
    if (this->p_alignreader != NULL){
        p_alignreader = this->p_alignreader;
        //throw runtime_error("p_alignreader == NULL");
    }
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
    vector<bool> var_data_temp(max_code + 4, false);
    for (int64_t i = 0; i < var_data.size(); ++i)
        var_data_temp[var_data[i].code] = true;
    
    
    // scan m5_file and recode
    if (p_alignreader==NULL)
        throw runtime_error("AlignCoderSNV::recode(): p_alignreader has not be set.");
    ofstream p_outfile;
    open_outfile(p_outfile, recode_file);
    ofstream p_outfile_ref;
    if (is_report_ref)
        open_outfile(p_outfile_ref, recode_file + ".ref");
    
    p_alignreader->open(m5_file);
    Align align;
    int nline = 0;
    while(p_alignreader->readline(align)){
        ++nline;
        if (nline % 100 == 0)
            cout << nline << endl;
        //cout << nline << endl;
        
        // exceptions
        int alen = (int) align.matchPattern.size();
        if ( !(align.qAlignedSeq.size()==alen && align.tAlignedSeq.size()==alen) )
            throw runtime_error("incorrect match patter in line " + to_string(nline));
        if (align.qStrand != '+')
            throw runtime_error("qStrand should be + in line " + to_string(nline));
        
        // reverse alignment if it is aligned to negative strand
        /*if (align.tStrand != '+'){
            align.qAlignedSeq = getrevcomp(align.qAlignedSeq);
            align.tAlignedSeq = getrevcomp(align.tAlignedSeq);
        }*/
        
        // encode
        int cur_pos = align.tStart;
        for (int i=0; i<alen; i++){
            //cout << "nline=" << nline << ", i=" <<i << endl;
            if (align.tAlignedSeq[i]=='-')
                continue;
            
            /*if (!(align.qAlignedSeq[i]=='A' || align.qAlignedSeq[i]=='C' || align.qAlignedSeq[i]=='G' || align.qAlignedSeq[i]=='T'))
                continue;*/
            
            if (4*cur_pos+3 > max_code + 3)
                break;
            
            // realign if hit detected variants
            int score_A = MIN_SCORE;
            int score_C = MIN_SCORE;
            int score_G = MIN_SCORE;
            int score_T = MIN_SCORE;
            bool is_var = false;

            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_A;
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_C;
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_G;
            seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> cur_realign_T;
            
            string cur_qseq;
            string cur_rseq;
            pair<string, string> context;
            char ref_base;
            
            // align local sequence to the referece
            if (var_data_temp[4*cur_pos] || var_data_temp[4*cur_pos+1] || var_data_temp[4*cur_pos+2] || var_data_temp[4*cur_pos+3]){
                //bool rl = this->get_context_m5(i, left_len, right_len, align.tAlignedSeq, context);
                bool rl = this->get_context_m5(i, left_len, right_len, align.tAlignedSeq, context, true);
                if (!rl){
                    ++cur_pos;
                    continue;
                }
                
                is_var = true;

                // get left query sequence length
                int64_t k = 0;
                int64_t cur_qseq_start = i;
                while(true){
                    if (align.tAlignedSeq[cur_qseq_start]!='-')
                        k++;
                    if (k >= context.first.size())
                        break;
                    --cur_qseq_start;
                }
                
                
                // get right query sequence length
                k = 0;
                int64_t cur_qseq_end = i+1;
                while(true){
                    if (align.tAlignedSeq[cur_qseq_end]!='-')
                        k++;
                    if (k >= context.second.size())
                        break;
                    ++cur_qseq_end;
                }
                
                if (cur_qseq_start < 0){
                    cur_qseq_start = 0;
                    //throw runtime_error("cur_qseq_start < 0");
                }
                if (cur_qseq_end >= alen){
                    cur_qseq_end = alen - 1;
                    //throw runtime_error("cur_qseq_end >= alen");
                }
                
                for (auto j = cur_qseq_start; j <= cur_qseq_end; ++j){
                    if (align.qAlignedSeq[j]!='-')
                        cur_qseq.push_back(align.qAlignedSeq[j]);
                }
                
                cur_rseq = context.first + context.second;
                
                if (cur_qseq == ""){
                    ++cur_pos;
                    continue;
                }
                
                if (cur_rseq == "")
                    throw runtime_error("cur_rseq is empty");
                
                ref_base = cur_rseq[context.first.size()-1];
                
                // realign
                cur_rseq[context.first.size()-1] = 'A';
                score_A = this->realign(cur_realign_A, cur_qseq, cur_rseq);
                
                cur_rseq[context.first.size()-1] = 'C';
                score_C = this->realign(cur_realign_C, cur_qseq, cur_rseq);
                
                cur_rseq[context.first.size()-1] = 'G';
                score_G = this->realign(cur_realign_G, cur_qseq, cur_rseq);
                
                cur_rseq[context.first.size()-1] = 'T';
                score_T = this->realign(cur_realign_T, cur_qseq, cur_rseq);
                
                // to be removed
                /*if (cur_pos==229419 && nline == 4){
                    cout << "cur_pos = " << cur_pos << endl;
                    cout << "A: " << score_A << endl << cur_realign_A;
                    cout << "C: " << score_C << endl << cur_realign_C;
                    cout << "G: " << score_G << endl << cur_realign_G;
                    cout << "T: " << score_T << endl << cur_realign_T;
                    int tmp = 0;
                }*/
            }else{
                ++cur_pos;
                continue;
            }
            
            // recode
            if (score_A == MIN_SCORE && score_C == MIN_SCORE && score_G == MIN_SCORE && score_T == MIN_SCORE)
                throw runtime_error("A,C,G,T == MIN_SCORE, no alignment was done");
            
            bool is_legal = true;
            //if (!(align.qAlignedSeq[i]=='A' || align.qAlignedSeq[i]=='C' || align.qAlignedSeq[i]=='G' || align.qAlignedSeq[i]=='T'))
            //    is_legal = false;
            // A
            if (score_A > score_C && score_A > score_G && score_A > score_T && is_legal){
                if (align.tAlignedSeq[i]!='A'){
                    p_outfile << 4*cur_pos << '\t';
                }else{
                    if (is_report_ref)
                        p_outfile_ref << 4*cur_pos << '\t';
                }
            }
            
            // C
            if (score_C > score_A && score_C > score_G && score_C > score_T && is_legal){
                if (align.tAlignedSeq[i]!='C'){
                    p_outfile << 4*cur_pos+1 << '\t';
                }else{
                    if (is_report_ref)
                        p_outfile_ref << 4*cur_pos+1 << '\t';
                }
            }
            
            // G
            if (score_G > score_A && score_G > score_C && score_G > score_T && is_legal){
                if (align.tAlignedSeq[i]!='G'){
                    p_outfile << 4*cur_pos+2 << '\t';
                }else{
                    if (is_report_ref)
                        p_outfile_ref << 4*cur_pos+2 << '\t';
                }
            }
            
            // T
            if (score_T > score_A && score_T > score_C && score_T > score_G && is_legal){
                if (align.tAlignedSeq[i]!='T'){
                    p_outfile << 4*cur_pos+3 << '\t';
                }else{
                    if (is_report_ref)
                        p_outfile_ref << 4*cur_pos+3 << '\t';
                }
                
            }
            
            ++cur_pos;
        }
        p_outfile << endl;
        if (is_report_ref)
            p_outfile_ref << endl;
    }
    
    p_alignreader->close();
    
    
    p_outfile.close();
    if (is_report_ref)
        p_outfile_ref.close();
    
    return true;
}

bool AlignCoderSNV::recode_multithread(string m5_file, string var_file, string recode_file, int left_len, int right_len, int nthread, bool is_report_ref)
{
    // split read files
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, m5_file);
    int64_t nlines = ceil ( (double)reads_range.size() / nthread);
    split_file(m5_file, recode_file + ".tmp.split", nlines);
    
    // multithread recode
    vector<thread> threads;
    for (auto i = 0; i < nthread; ++i){
        cout << "recode thread " << i + 1 << endl;
        string cur_m5file = recode_file + ".tmp.split.part." + to_string(i);
        string cur_recodefile = recode_file + ".tmp.split.part." + to_string(i) + ".recode";
        cout << cur_m5file << endl;
        cout << cur_recodefile << endl;
        threads.push_back(thread(&AlignCoderSNV::recode, this, cur_m5file, var_file, cur_recodefile, left_len, right_len, is_report_ref));
    }
    
    for (auto i = 0; i < threads.size(); ++i)
        threads[i].join();
    
    // merge results
    string cmd_merge_recode = "cat ";
    string cmd_merge_recode_ref = "cat ";
    for (auto i = 0; i < nthread; ++i){
        cmd_merge_recode += recode_file + ".tmp.split.part." + to_string(i) + ".recode ";
        cmd_merge_recode_ref += recode_file + ".tmp.split.part." + to_string(i) + ".recode.ref ";
    }
    cmd_merge_recode += "> " + recode_file;
    cmd_merge_recode_ref += "> " + recode_file + ".ref";
    
    cout << cmd_merge_recode << endl; system(cmd_merge_recode.c_str());
    cout << cmd_merge_recode_ref << endl; system(cmd_merge_recode_ref.c_str());
    
    // remove temp results
    string cmd = "rm " + recode_file + ".tmp.split.part.*";
    cout << cmd << endl; system(cmd.c_str());
    
    return true;
}

bool AlignCoderSNV::get_context_m5(int i, int left, int right, const string &tAlignedSeq, pair<string,string> &context, bool is_overhanged)
{
    int n_base, k;
    
    // search left
    n_base = 0; k = 0;
    char cur_base = tAlignedSeq[i];
    string context_left("");
    while (n_base<=left){
        k++;
        if (i - k < 0){
            if (is_overhanged)
                break;
            else
                return false;
        }
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
        if (i + k >= tAlignedSeq.size()){
            if (is_overhanged)
                break;
            else
                return false;
        }
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


int AlignCoderSNV::realign(seqan::Align<seqan::String<seqan::Dna5>, seqan::ArrayGaps> &cur_realign, string _qseq, string _rseq)
{
    seqan::String<seqan::Dna5> qseq = _qseq;
    seqan::String<seqan::Dna5> rseq = _rseq;
    
    seqan::resize(seqan::rows(cur_realign),2);
    seqan::assignSource(seqan::row(cur_realign, 0), qseq);
    seqan::assignSource(seqan::row(cur_realign, 1), rseq);
    
    int gapOpenScore = -4;
    int gapExtendScore = -2;
    
    seqan::Score<int, seqan::ScoreMatrix<seqan::Dna5, seqan::Default> > cur_score(gapExtendScore, gapOpenScore);
    
    for (auto i = 0; i < 4; ++i)
    {
        for (auto j = 0; j < 4; ++j)
        {
            if (i == j)
                setScore(cur_score, seqan::Dna5(i), seqan::Dna5(j), 2);
            else
                setScore(cur_score, seqan::Dna5(i), seqan::Dna5(j), -4);
        }
    }
    for (auto i = 0; i < 5; ++i){
        setScore(cur_score, seqan::Dna5(i), seqan::Dna5(4), 0);
        setScore(cur_score, seqan::Dna5(4), seqan::Dna5(i), 0);
    }
    
    int cur_align_score = seqan::globalAlignment(cur_realign, cur_score, seqan::AffineGaps());
    
    //int cur_align_score = seqan::globalAlignment(cur_realign, seqan::Score<int, seqan::Simple>(2, -4, -2, -4), seqan::AffineGaps());
    return cur_align_score;
}








