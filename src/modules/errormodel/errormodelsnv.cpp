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
    
    // get context effect without summary
    cout << "get context effect without summary" << endl;
    ContextEffectAll context_effect_all;
    get_context_effect_all(pileup, context_effect_all);
    print_context_effect_all(out_prefix + ".context.all", context_effect_all);
    
    // get context effect summary
    cout << "get context effect" << endl;
    ContextEffect context_effect;
    get_context_effect(pileup, context_effect);
    print_context_effect(out_prefix + ".context", context_effect);
    
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
        if (align.tStrand != '+')
            throw runtime_error("there are reads aligned to negative strand.");
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
                    if (align.tAlignedSeq[i] != pileup[locus].ref && pileup[locus].ref!='$'){
                        cout << "reads: " << align.qName << "," << i << endl;
                        cout << "locus: " << locus << endl;
                        cout << "align.tAlignedSeq[i]: " << align.tAlignedSeq[i] << ", ref: " << pileup[locus].ref << endl;
                        print_pileup("../results/B_10_cons.pileup", pileup);
                        throw runtime_error("align.tAlignedSeq[i] != pileup[locus].ref");
                    }
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
                // get ref
                pileup[locus].ref = align.tAlignedSeq[i];
            }
            locus++;
        }
        n_reads++;
        //break;
        if (n_reads % 1000 == 0)
            cout << "processed " << n_reads << " reads" << endl;
    }
    cout << "processed " << n_reads << " reads" << endl;
    alignreader.close();
}

void ErrorModelSNV::print_pileup(string out_file, const vector<BaseFreq> &pileup){
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)pileup.size(); i++){
        fs_out << i+1 << "\t" << pileup[i].context.first << "," << pileup[i].context.second << '\t';
        fs_out << pileup[i].ref << '\t';
        fs_out << "A:" << pileup[i].nvar[0] << ',';
        fs_out << "C:" << pileup[i].nvar[1] << ',';
        fs_out << "G:" << pileup[i].nvar[2] << ',';
        fs_out << "T:" << pileup[i].nvar[3] << '\t';
        fs_out << pileup[i].cvg << endl;
    }
    fs_out.close();
}

void ErrorModelSNV::get_context_effect_all(const vector<BaseFreq> &pileup, ContextEffectAll &context_effect_all)
{
    for (int i=0; i<(int) pileup.size(); i++){
        if (pileup[i].context.first=="" || pileup[i].context.second=="")
            continue;
        context_effect_all[pileup[i].context.first][pileup[i].context.second].push_back(pileup[i]);
    }
}

void ErrorModelSNV::get_context_effect(const vector<BaseFreq> &pileup, ContextEffect &context_effect)
{
    for (int i=0; i<(int) pileup.size(); i++){
        if (pileup[i].context.first=="" || pileup[i].context.second=="")
            continue;
        context_effect[pileup[i].context.first][pileup[i].context.second].cvg += pileup[i].cvg;
        context_effect[pileup[i].context.first][pileup[i].context.second].nvar[0] += pileup[i].nvar[0];
        context_effect[pileup[i].context.first][pileup[i].context.second].nvar[1] += pileup[i].nvar[1];
        context_effect[pileup[i].context.first][pileup[i].context.second].nvar[2] += pileup[i].nvar[2];
        context_effect[pileup[i].context.first][pileup[i].context.second].nvar[3] += pileup[i].nvar[3];
    }
}

void ErrorModelSNV::print_context_effect_all(string out_file, ContextEffectAll &context_effect_all)
{
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    unordered_map<string, unordered_map<string, vector<BaseFreq> > >::iterator it_1;
    unordered_map<string, vector<BaseFreq> >::iterator it_2;
    for (it_1=context_effect_all.begin(); it_1!=context_effect_all.end(); it_1++){
        for (it_2 = it_1->second.begin(); it_2!=it_1->second.end(); it_2++){
            fs_out << it_1->first << "," << it_2->first << '\t';
            for (int i=0; i<(int)it_2->second.size(); i++){
                fs_out << "A:" << it_2->second[i].nvar[0] << ',';
                fs_out << "C:" << it_2->second[i].nvar[1] << ',';
                fs_out << "G:" << it_2->second[i].nvar[2] << ',';
                fs_out << "T:" << it_2->second[i].nvar[3] << ',';
                fs_out << "cvg:" << it_2->second[i].cvg << '\t';
            }
            fs_out << endl;
        }
    }
    fs_out.close();
}
void ErrorModelSNV::print_context_effect(string out_file, ContextEffect &context_effect)
{
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    unordered_map<string, unordered_map<string, BaseFreq > >::iterator it_1;
    unordered_map<string, BaseFreq >::iterator it_2;
    for (it_1=context_effect.begin(); it_1!=context_effect.end(); it_1++){
        for (it_2 = it_1->second.begin(); it_2!=it_1->second.end(); it_2++){
            fs_out << it_1->first << '\t' << it_2->first << '\t';
            fs_out << it_2->second.nvar[0] << '\t';
            fs_out << it_2->second.nvar[1] << '\t';
            fs_out << it_2->second.nvar[2] << '\t';
            fs_out << it_2->second.nvar[3] << '\t';
            fs_out << it_2->second.cvg << endl;
        }
    }
    fs_out.close();
}

// only work for forward reads
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


