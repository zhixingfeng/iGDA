//
//  errormodelsnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/10.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "errormodelsnv.h"

void ErrorModelSNV::learn(string align_file, string out_prefix)
{
    // get size of genome covered by reads
    cout << "get genome size" << endl;
    int g_size = get_genomesize(align_file);
    vector<BaseFreq> pileup(g_size, BaseFreq());
    cout << "genome size is " << g_size << endl;
    
    // pileup 
    cout << "pileup reads" << endl;
    pileup_reads(align_file, pileup);
    print_pileup(out_prefix + ".pileup", pileup);
    
    // get context effect
    
}

int ErrorModelSNV::get_genomesize(string align_file){
    // tEnd in align is 0-based, should add 1 in the end!
    int g_size = -1;
    AlignReaderM5 alignreader;
    alignreader.open(align_file);
    while (1){
        Align align;
        if (!alignreader.readline(align))
            break;
        if (align.tEnd > g_size)
            g_size = align.tEnd;
    }
    alignreader.close();
    return g_size + 1;
}


void ErrorModelSNV::pileup_reads(string align_file, vector<BaseFreq> &pileup)
{
    
    AlignReaderM5 alignreader;
    alignreader.open(align_file);
    int64_t n_reads = 0;
    while (1){
        Align align;
        if (!alignreader.readline(align))
            break;
        if (align.mapQV!=254)
            continue;
        int locus = align.tStart;    
        for (int i=0; i<(int)align.tAlignedSeq.size(); i++){
            // get context
            pair<string, string> context;
            if (align.tAlignedSeq[i]=='-')
                continue;
            
            if (get_context_m5(i, 1, 1, align.tAlignedSeq, context)){
                pileup[locus].context = context;

                // add coverage 
                pileup[locus].cvg++;
            
                // get SNV
                if (align.qAlignedSeq[i] != align.tAlignedSeq[i] && align.qAlignedSeq[i]!='-'){
                    switch (align.qAlignedSeq[i]) {
                        case 'A':
                            pileup[locus].nvar[0]++;
                            break;
                        case 'C':
                            pileup[locus].nvar[1]++;
                            break;
                        case 'G':
                            pileup[locus].nvar[2]++;
                            break;
                        case 'T':
                            pileup[locus].nvar[3]++;
                            break;
                        default:
                            throw runtime_error("sequence should not contain bases other than A, C, G, T");
                            break;
                    }
                }
            }
            locus++;
        }
        n_reads++;
        if (n_reads % 1000 == 0)
            cout << "processed " << n_reads << " reads" << endl;
    }
    alignreader.close();
}

void ErrorModelSNV::print_pileup(string out_file, const vector<BaseFreq> &pileup){
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)pileup.size(); i++){
        fs_out << i+1 << "\t" << pileup[i].context.first << "," << pileup[i].context.second << '\t';
        fs_out << "A:" << pileup[i].nvar[0] << ',';
        fs_out << "C:" << pileup[i].nvar[1] << ',';
        fs_out << "G:" << pileup[i].nvar[2] << ',';
        fs_out << "T:" << pileup[i].nvar[3] << '\t';
        fs_out << pileup[i].cvg << endl;
    }
    fs_out.close();
}

bool ErrorModelSNV::get_context_m5(int i, int left, int right, const string &tAlignedSeq, pair<string,string> &context)
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


