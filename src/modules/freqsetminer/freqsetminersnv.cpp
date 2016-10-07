//
//  freqsetminersnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 9/26/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "freqsetminersnv.h"


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
