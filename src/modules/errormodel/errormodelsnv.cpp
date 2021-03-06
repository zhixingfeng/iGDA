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

void ErrorModelSNV::merge(vector<string> &context_files)
{
    // merge context effect
    ContextEffect context_effect;
    for (int i=0; i<(int)context_files.size(); i++){
        ifstream fs_in;
        open_infile(fs_in, context_files[i]);
        while (1){
            BaseFreq buf;
            fs_in >> buf.context.first >> buf.context.second;
            fs_in >> buf.nvar[0] >> buf.nvar[1] >> buf.nvar[2] >> buf.nvar[3];
            fs_in >> buf.cvg;
            if (fs_in.eof())
                break;
            context_effect[buf.context.first][buf.context.second].cvg += buf.cvg;
            context_effect[buf.context.first][buf.context.second].nvar[0] += buf.nvar[0];
            context_effect[buf.context.first][buf.context.second].nvar[1] += buf.nvar[1];
            context_effect[buf.context.first][buf.context.second].nvar[2] += buf.nvar[2];
            context_effect[buf.context.first][buf.context.second].nvar[3] += buf.nvar[3];
        }
        fs_in.close();
    }
    
    // print context effect
    unordered_map<string, unordered_map<string, BaseFreq > >::iterator it_1;
    unordered_map<string, BaseFreq >::iterator it_2;
    for (it_1=context_effect.begin(); it_1!=context_effect.end(); it_1++){
        for (it_2 = it_1->second.begin(); it_2!=it_1->second.end(); it_2++){
            cout << it_1->first << '\t' << it_2->first << '\t';
            cout << it_2->second.nvar[0] << '\t';
            cout << it_2->second.nvar[1] << '\t';
            cout << it_2->second.nvar[2] << '\t';
            cout << it_2->second.nvar[3] << '\t';
            cout << it_2->second.cvg << endl;
        }
    }

}

void ErrorModelSNV::merge_all(vector<string> &context_all_files)
{
    // merge context effect all 
    ContextEffectAll context_effect_all;
    for (int i=0; i<(int)context_all_files.size(); i++){
        ifstream fs_infile;
        open_infile(fs_infile, context_all_files[i]);

        while (1){
            if (!readline_context_effect_all(context_effect_all, fs_infile))
                break;
        }
        
        fs_infile.close();
    }
    
    // print context effect all 
    unordered_map<string, unordered_map<string, vector<BaseFreq> > >::iterator it_1;
    unordered_map<string, vector<BaseFreq> >::iterator it_2;
    for (it_1=context_effect_all.begin(); it_1!=context_effect_all.end(); it_1++){
        for (it_2 = it_1->second.begin(); it_2!=it_1->second.end(); it_2++){
            cout << it_1->first << "," << it_2->first << '\t';
            for (int i=0; i<(int)it_2->second.size(); i++){
                cout << "A:" << it_2->second[i].nvar[0] << ',';
                cout << "C:" << it_2->second[i].nvar[1] << ',';
                cout << "G:" << it_2->second[i].nvar[2] << ',';
                cout << "T:" << it_2->second[i].nvar[3] << ',';
                cout << "cvg:" << it_2->second[i].cvg << '\t';
            }
            cout << endl;
        }
    }

}



void ErrorModelSNV::pileup_count_to_context(string pu_count_file, string pu_file, string context_file)
{
    // load context for each locus
    unordered_map<int64_t, pair<string, string> > locus_context;
    unordered_map<int64_t, char > locus_ref;
    ifstream fs_pu_file;
    open_infile(fs_pu_file, pu_file);
    while(true){
        string buf;
        getline(fs_pu_file, buf);
        if(fs_pu_file.eof())
            break;
        
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size() != 9)
            throw runtime_error("incorrect format in " + pu_file);
        
        int64_t locus = int64_t(stod(buf_vec[0]));
        string context_left = buf_vec[1];
        string context_right = buf_vec[2];
        char ref = buf_vec[3][0];
        locus_context[locus] = pair<string, string>(context_left, context_right);
        locus_ref[locus] = ref;
        
    }
    fs_pu_file.close();
    
    // get context effect
    ContextEffect context_effect;
    ifstream fs_pu_count_file;
    open_infile(fs_pu_count_file, pu_count_file);
    int64_t prev_locus = -1;
    while(true){
        // read line
        string buf;
        getline(fs_pu_count_file, buf);
        if(fs_pu_count_file.eof())
            break;
        
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size() != 6)
            throw runtime_error("incorrect format in " + pu_file);
        
        // get context effect
        int64_t locus = int64_t(stod(buf_vec[0]));
        char var_base = buf_vec[1][0];
        int64_t n_var = int64_t(stod(buf_vec[3]));
        int64_t cvg = int64_t(stod(buf_vec[4]));
    
        auto it = locus_context.find(locus);
        if (it==locus_context.end())
            throw runtime_error("fail to find context");
        
        string context_left = it->second.first;
        string context_right = it->second.second;
        
        auto it_2 = locus_ref.find(locus);
        if (it_2 == locus_ref.end())
            throw runtime_error("fail to find ref");
        
        char cur_ref = it_2->second;
        
        // fill context effect
        auto it_left = context_effect.find(context_left);
        if (it_left == context_effect.end()){
            context_effect[context_left][context_right].cvg = 0;
            context_effect[context_left][context_right].nvar[0] = 0;
            context_effect[context_left][context_right].nvar[1] = 0;
            context_effect[context_left][context_right].nvar[2] = 0;
            context_effect[context_left][context_right].nvar[3] = 0;
        }else{
            auto it_right = it_left->second.find(context_right);
            if (it_right == it_left->second.end()){
                context_effect[context_left][context_right].cvg = 0;
                context_effect[context_left][context_right].nvar[0] = 0;
                context_effect[context_left][context_right].nvar[1] = 0;
                context_effect[context_left][context_right].nvar[2] = 0;
                context_effect[context_left][context_right].nvar[3] = 0;
            }
        }
        context_effect[context_left][context_right].ref = cur_ref;
        if (locus != prev_locus){
            context_effect[context_left][context_right].cvg += cvg;
            prev_locus = locus;
        }
        switch (var_base) {
            case 'A':
                context_effect[context_left][context_right].nvar[0] += n_var;
                break;
            
            case 'C':
                context_effect[context_left][context_right].nvar[1] += n_var;
                break;
            
            case 'G':
                context_effect[context_left][context_right].nvar[2] += n_var;
                break;
            
            case 'T':
                context_effect[context_left][context_right].nvar[3] += n_var;
                break;
            
            default:
                throw runtime_error("incorrect var_base");
                break;
        }
        
        
    }
    
    fs_pu_count_file.close();
    
    this->print_context_effect(context_file, context_effect);
}







int ErrorModelSNV::get_genomesize(string align_file)
{
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
    bool is_neg_strand = false;
    int64_t n_reads = 0;
    while (1){
        Align align;
        if (!alignreader.readline(align))
            break;
        if (align.tStrand != '+')
            is_neg_strand = true;
        //    throw runtime_error("there are reads aligned to negative strand.");
        //if (align.mapQV!=254)
        //    continue;
        
        int locus = align.tStart;    
        for (int i=0; i<(int)align.tAlignedSeq.size(); i++){
            // get context
            pair<string, string> context;
            if (align.tAlignedSeq[i]=='-')
                continue;
            
            if (get_context_m5(i, this->left_len, this->right_len, align.tAlignedSeq, context)){
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
                        case 'N':
                            break;
                        default:
                            throw runtime_error("sequence should not contain bases other than A, C, G, T or N");
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
    
    if (is_neg_strand)
        cout << "Warning: input file contains reads aligned to negative strand";
}

void ErrorModelSNV::print_pileup(string out_file, const vector<BaseFreq> &pileup){
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)pileup.size(); i++){
        if (pileup[i].context.first=="" || pileup[i].context.second=="" || pileup[i].ref=='$')
            continue;
        fs_out << i << '\t' << pileup[i].context.first << '\t' << pileup[i].context.second << '\t';
        fs_out << pileup[i].ref << '\t';
        fs_out << pileup[i].nvar[0] << '\t';  // A
        fs_out << pileup[i].nvar[1] << '\t';  // C
        fs_out << pileup[i].nvar[2] << '\t';  // G
        fs_out << pileup[i].nvar[3] << '\t';  // T
        fs_out << pileup[i].cvg << endl;
    }
    fs_out.close();
}

vector<BaseFreq> ErrorModelSNV::load_pileup(string pu_file)
{
    // scan the pileup file to get the genome_size;
    int g_size = -1;
    ifstream fs_in;
    open_infile(fs_in, pu_file);
    while (true) {
        string buf;
        getline(fs_in, buf);
        
        if (fs_in.eof()) break;
        
        vector<string> buf_vec = split(buf, '\t');
        
        if (buf_vec.size() != 9)
            throw runtime_error("load_pileup(): buf_vec.size() != 9");
        
        int locus = stoi(buf_vec[0]);
        if (locus > g_size)
            g_size = locus;
    }
    ++g_size;
    fs_in.close();
    
    // load pileup
    vector<BaseFreq> pileup(g_size, BaseFreq());
    open_infile(fs_in, pu_file);
    while (true) {
        string buf;
        getline(fs_in, buf);
        
        if (fs_in.eof()) break;
        
        vector<string> buf_vec = split(buf, '\t');
        
        if (buf_vec.size() != 9)
            throw runtime_error("load_pileup(): buf_vec.size() != 9");
        
        if (buf_vec[3].size() != 1)
            throw runtime_error("load_pileup(): buf_vec[3].size() != 1");
        
        int locus = stoi(buf_vec[0]);
        BaseFreq basefreq;
        basefreq.context.first = buf_vec[1];
        basefreq.context.second = buf_vec[2];
        basefreq.ref = buf_vec[3][0];
        basefreq.nvar[0] = stoi(buf_vec[4]);
        basefreq.nvar[1] = stoi(buf_vec[5]);
        basefreq.nvar[2] = stoi(buf_vec[6]);
        basefreq.nvar[3] = stoi(buf_vec[7]);
        basefreq.cvg = stoi(buf_vec[8]);
        pileup[locus] = basefreq;
    }
    //cout << g_size << endl;
    return pileup;
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

bool ErrorModelSNV::readline_context_effect_all(ContextEffectAll &context_effect_all, ifstream & fs_infile)
{
    string buf;
    getline(fs_infile, buf);
    if (fs_infile.eof())
        return false;
    vector<string> buf_split = split(buf, '\t');
    if (buf_split.size()<2)
        throw runtime_error("readline_context_effect_all: incorrect format");
    
    // get context
    vector<string> cur_context = split(buf_split[0], ',');
    if (cur_context.size()!=2)
        throw runtime_error("incorrect cur_context");

    // add positions with the same context
    for (int i=1; i<(int)buf_split.size(); i++){
        vector<string> cur_freq_list = split(buf_split[i], ',');
        if (cur_freq_list.size()!=5)
            throw runtime_error("incorrect format cur_freq_list");
        BaseFreq cur_freq;
        vector<string> tmp;
        
        // load freq of A
        tmp = split(cur_freq_list[0],':');
        if (tmp.size()!=2)
            throw runtime_error("incorrect format cur_freq_list[0]");
        cur_freq.nvar[0] = stoi(tmp[1]);
        
        // load freq of C
        tmp = split(cur_freq_list[1],':');
        if (tmp.size()!=2)
            throw runtime_error("incorrect format cur_freq_list[1]");
        cur_freq.nvar[1] = stoi(tmp[1]);
        
        // load freq of G
        tmp = split(cur_freq_list[2],':');
        if (tmp.size()!=2)
            throw runtime_error("incorrect format cur_freq_list[2]");
        cur_freq.nvar[2] = stoi(tmp[1]);
        
        // load freq of T
        tmp = split(cur_freq_list[3],':');
        if (tmp.size()!=2)
            throw runtime_error("incorrect format cur_freq_list[3]");
        cur_freq.nvar[3] = stoi(tmp[1]);
        
        // load cvg 
        tmp = split(cur_freq_list[4],':');
        if (tmp.size()!=2)
            throw runtime_error("incorrect format cur_freq_list[4]");
        cur_freq.cvg = stoi(tmp[1]);

        // add cur_freq to context_effect_all
        context_effect_all[cur_context[0]][cur_context[1]].push_back(cur_freq);
    }
    
    
    return true;
}


