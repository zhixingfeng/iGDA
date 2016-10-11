//
//  freqsetminersnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 9/26/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "freqsetminersnv.h"

vector<int> FreqSetMinerSNV::detect(string encode_file, string align_file, string cmpreads_file, vector<double> p0, double log_bf_cutoff)
{
    vector<int> var_loci;
    
    map<int,int> var_loci_map;
    //map<int,int>::iterator it;
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    // pile up
    vector<int> pu = pileup(encode_file);
    
    
    ifstream fs_cmpreads_file; open_infile(fs_cmpreads_file, cmpreads_file);
    int n=1;
    while (1){
        if (n%100==0)
            cout << n << endl;
        n++;
        CmpReads cmpreads_line;
        if (!readLineCmpReads(fs_cmpreads_file, cmpreads_line))
            break;
    
        vector<CondFreq> condfreq = getcondfreq(cmpreads_line.cons_code, encode_data, reads_range, pu, p0);
        if (condfreq.size()!=cmpreads_line.cons_code.size())
            throw runtime_error("condfreq.size()!=cmpreads_line.cons_code.size()");
        for (int i=0; i<(int)condfreq.size(); i++){
            if (condfreq[i].log_bf>=log_bf_cutoff){
                var_loci_map[cmpreads_line.cons_code[i]]++;
            }
        }
    }
    cout << n << endl;
    fs_cmpreads_file.close();
    
    map<int,int>::iterator it;
    for (it=var_loci_map.begin(); it!=var_loci_map.end(); it++)
        var_loci.push_back(it->first);
    
    return var_loci;
}

vector<int> FreqSetMinerSNV::detectVariantsCoarse(string encode_file, string align_file, string cmpreads_file, double p_cutoff)
{
    // set mutation rate p_mut, should be replaced by better estimation
    double p_mut = 0.01;
    
    if (ptr_aligncoder==NULL)
        throw runtime_error("ptr_aligncoder is NULL");

    // pileup and caclulate coverage
    vector<int> pu = pileup(encode_file);
    vector<int> cvg = getcvg(align_file);
    if (ptr_aligncoder->decode((int)pu.size()).first > cvg.size())
        throw runtime_error("pu is larger than cvg");
    
    vector<int> rl;
    ifstream fs_cmpreads_file; open_infile(fs_cmpreads_file, cmpreads_file);
    
    vector<double> p_values;
    int n=1;
    while(1){
        if (n%10000==0)
            cout << n << endl;
        CmpReads cmpreads_line;
        if (readLineCmpReads(fs_cmpreads_file, cmpreads_line)==false)
            break;
        
        double stat = 0;
        double cvg_sum = 0;
        for (int i=0;i<(int)cmpreads_line.cons_code.size(); i++){
            int cur_locus =ptr_aligncoder->decode(cmpreads_line.cons_code[i]).first;
            if (cmpreads_line.cons_code[i]>pu.size())
                throw runtime_error("cmpreads_line.cons_code[i]>pu.size()");
            
            if (cur_locus>cvg.size())
                throw runtime_error("cur_locus>cvg.size()");

            stat += pu[cmpreads_line.cons_code[i]-1];
            cvg_sum += cvg[cur_locus-1];
            
        }
        double mu = cvg_sum*p_mut;
        double sigma = sqrt(mu*(1-p_mut));
        double z_score = (stat - mu)/sigma;
        p_values.push_back(1 - pnorm(z_score));
        
        n++;
    }
    cout << n << endl;
    
    fs_cmpreads_file.close();
    return rl;
}

vector<double> FreqSetMinerSNV::detectVariantsSingle(string encode_file, string align_file)
{
    // set mutation rate p_mut, should be replaced by better estimation
    double p_mut = 0.01;
    
    if (ptr_aligncoder==NULL)
        throw runtime_error("ptr_aligncoder is NULL");
    vector<double> p_values;
    vector<int> pu = pileup(encode_file);
    vector<int> cvg = getcvg(align_file);
    
    if (ptr_aligncoder->decode((int)pu.size()).first > cvg.size())
        throw runtime_error("pu is larger than cvg");
    
    // scan the genome and test
    for (int i=0; i<(int)pu.size(); i++){
        int locus = ptr_aligncoder->decode(i+1).first;
        int cur_cvg = cvg[locus-1];
        
        double mu = cur_cvg * p_mut;
        double sigma = sqrt(mu * (1 - p_mut));
        double z_score = (pu[i] - mu) / sigma;
        p_values.push_back( 1 - pnorm(z_score) );
     }
    
    return p_values;
}

void FreqSetMinerSNV::getMarginalFreq(string encode_file, string align_file, string out_file)
{
    if (ptr_aligncoder==NULL)
        throw runtime_error("ptr_aligncoder is NULL");

    vector<int> pu = pileup(encode_file);
    vector<int> cvg = getcvg(align_file);
    
    if (ptr_aligncoder->decode((int)pu.size()).first > cvg.size())
        throw runtime_error("pu is larger than cvg");
    
    // calculate frquency
    vector<double> prob(pu.size(),-1);
    vector<int> cvg_vec(pu.size(),-1);
    double avg_prob = 0;
    double n = 0;
    for (int i=0; i<(int)pu.size(); i++){
        int locus = ptr_aligncoder->decode(i+1).first;
        int cur_cvg = cvg[locus-1];
        cvg_vec[i] = cur_cvg;
        if (cur_cvg==0)
            prob[i] = -1;
        else{
            prob[i] = 1.0*pu[i] / cur_cvg;
            avg_prob += prob[i];
            n++;
        }
    }
    
    // fill 0 and -1 with avg_prob
    if (n==0) throw runtime_error("n==0");
    avg_prob = avg_prob / n;
    for (int i=0; i<(int)prob.size(); i++)
        if (prob[i]==0 || prob[i]==-1)
            prob[i] = avg_prob;
    
    // print prob
    ofstream fs_out_file; open_outfile(fs_out_file, out_file);
    for (int i=0; i<(int)prob.size(); i++)
        fs_out_file << prob[i] << '\t' << pu[i] << '\t' << cvg_vec[i] << endl;
    fs_out_file.close();
}


