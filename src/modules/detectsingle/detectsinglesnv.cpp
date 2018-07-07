//
//  detectsinglesnv.cpp
//  iGDA
//
//  Created by Zhixing Feng on 7/5/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "detectsinglesnv.h"


void DetectSingleSNV::loadcontexteffect(string contexteffect_file, int min_context_cvg)
{
    ifstream fs_infile;
    open_infile(fs_infile, contexteffect_file);
    while(true){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size() != 7)
            throw runtime_error("incorrect format in " + contexteffect_file);
        
        double A_freq = stod(buf_vec[2]);
        double C_freq = stod(buf_vec[3]);
        double G_freq = stod(buf_vec[4]);
        double T_freq = stod(buf_vec[5]);
        double cvg = stod(buf_vec[6]);
        
        if (cvg < min_context_cvg) continue;

        vector<double> &cur_context = this->contexteffect[buf_vec[0]][buf_vec[1]];
        if (cur_context.size() != 0)
            throw runtime_error("redundant context effect in " + contexteffect_file);
        
        cur_context.push_back(A_freq / cvg);
        cur_context.push_back(C_freq / cvg);
        cur_context.push_back(G_freq / cvg);
        cur_context.push_back(T_freq / cvg);
        cur_context.push_back(cvg);
        
        
        //cout << cur_context << endl;
    }
    
    fs_infile.close();
}

void DetectSingleSNV::savecontexteffect(string outfile)
{
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for (auto it_1 = this->contexteffect.begin(); it_1 != this->contexteffect.end(); ++it_1){
        for(auto it_2 = it_1->second.begin(); it_2 != it_1->second.end(); ++it_2){
            if (it_2->second.size() != 5)
                throw runtime_error("contexteffect should have 5 elements.");
            fs_outfile << it_1->first << '\t' << it_2->first << '\t';
            fs_outfile << int(it_2->second[0]*it_2->second[4]+0.5) << '\t';
            fs_outfile << int(it_2->second[1]*it_2->second[4]+0.5) << '\t';
            fs_outfile << int(it_2->second[2]*it_2->second[4]+0.5) << '\t';
            fs_outfile << int(it_2->second[3]*it_2->second[4]+0.5) << '\t';
            fs_outfile << int(it_2->second[4]+0.5) << endl;
        }
    }
    fs_outfile.close();
}

void DetectSingleSNV::detect(string pileup_file, string out_file, double min_bf, double min_prop, int min_cvg)
{
    if (this->contexteffect.size() == 0)
        throw runtime_error("DetectSingleSNV::detect: no context effect loaded.");
    
}




