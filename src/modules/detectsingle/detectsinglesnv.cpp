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
    
    ifstream fs_infile;
    open_infile(fs_infile, pileup_file);
    while(true){
        // load pileup record
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size() != 9)
            throw runtime_error("incorrect format in " + pileup_file);
        
        DetectSingleResult cur_result;
        cur_result.locus = stoi(buf_vec[0]);
        cur_result.context_left = buf_vec[1];
        cur_result.context_right = buf_vec[2];
        cur_result.ref_base = buf_vec[3][0];
        cur_result.A_count = stoi(buf_vec[4]);
        cur_result.C_count = stoi(buf_vec[5]);
        cur_result.G_count = stoi(buf_vec[6]);
        cur_result.T_count = stoi(buf_vec[7]);
        cur_result.cvg = stoi(buf_vec[8]);
        
        cur_result.A_freq = cur_result.cvg == 0 ? -1 : double(cur_result.A_count) / cur_result.cvg;
        cur_result.A_freq_ref = -1; cur_result.A_bf = -1;
        
        cur_result.C_freq = cur_result.cvg == 0 ? -1 :double(cur_result.C_count) / cur_result.cvg;
        cur_result.C_freq_ref = -1; cur_result.C_bf = -1;
        
        cur_result.G_freq = cur_result.cvg == 0 ? -1 :double(cur_result.G_count) / cur_result.cvg;
        cur_result.G_freq_ref = -1; cur_result.G_bf = -1;
        
        cur_result.T_freq = cur_result.cvg == 0 ? -1 :double(cur_result.T_count) / cur_result.cvg;
        cur_result.T_freq_ref = -1; cur_result.T_bf = -1;
        
        cur_result.cvg_ref = -1;
        
        if (cur_result.cvg == 0)
            continue;
        
        // calculate bayes factor based on context effect
        auto it_1 = this->contexteffect.find(cur_result.context_left);
        if (it_1 == this->contexteffect.end())
            continue;
        auto it_2 = it_1->second.find(cur_result.context_right);
        if (it_2 == it_1->second.end())
            continue;
        
        // test A
        if (it_2->second[0] > 0){
            cur_result.A_freq_ref = it_2->second[0];
            cur_result.A_bf = binom_log_bf(cur_result.A_count, cur_result.cvg, cur_result.A_freq_ref);
        }
        
        // test C
        if (it_2->second[1] > 0){
            cur_result.C_freq_ref = it_2->second[0];
            cur_result.C_bf = binom_log_bf(cur_result.C_count, cur_result.cvg, cur_result.C_freq_ref);
        }
        
        // test G
        if (it_2->second[2] > 0){
            cur_result.G_freq_ref = it_2->second[2];
            cur_result.G_bf = binom_log_bf(cur_result.G_count, cur_result.cvg, cur_result.G_freq_ref);
        }
        
        // test T
        if (it_2->second[3] > 0){
            cur_result.T_freq_ref = it_2->second[3];
            cur_result.G_bf = binom_log_bf(cur_result.G_count, cur_result.cvg, cur_result.G_freq_ref);
        }
        
        if (it_2->second[4] > 0)
            cur_result.cvg_ref = it_2->second[4];
        
        // store results
        this->result.push_back(cur_result);
    }
    
    fs_infile.close();
    
    // print results
    this->print_result(out_file);
}

void DetectSingleSNV::print_result(string out_file, double min_bf, double min_prop, int min_cvg)
{
    ofstream fs_outfile;
    open_outfile(fs_outfile, out_file);
    for (auto i = 0; i < this->result.size(); ++i){
        // check coverage
        if (this->result[i].cvg < min_cvg)
            continue;
        
        // check A
        if (this->result[i].A_freq >= min_prop && this->result[i].A_bf >= min_bf){
            fs_outfile << 4*this->result[i].locus << '\t';
            fs_outfile << this->result[i].A_bf << '\t';
            fs_outfile << this->result[i].A_freq << '\t';
            fs_outfile << this->result[i].A_count << '\t';
            fs_outfile << this->result[i].cvg << '\t';
            fs_outfile << 0 << '\t' << -1 << endl;
        }
        
        // check C
        if (this->result[i].C_freq >= min_prop && this->result[i].C_bf >= min_bf){
            fs_outfile << 4*this->result[i].locus << '\t';
            fs_outfile << this->result[i].C_bf << '\t';
            fs_outfile << this->result[i].C_freq << '\t';
            fs_outfile << this->result[i].C_count << '\t';
            fs_outfile << this->result[i].cvg << '\t';
            fs_outfile << 0 << '\t' << -1 << endl;
        }
        
        // check G
        if (this->result[i].G_freq >= min_prop && this->result[i].G_bf >= min_bf){
            fs_outfile << 4*this->result[i].locus << '\t';
            fs_outfile << this->result[i].G_bf << '\t';
            fs_outfile << this->result[i].G_freq << '\t';
            fs_outfile << this->result[i].G_count << '\t';
            fs_outfile << this->result[i].cvg << '\t';
            fs_outfile << 0 << '\t' << -1 << endl;
        }
        
        // check T
        if (this->result[i].T_freq >= min_prop && this->result[i].T_bf >= min_bf){
            fs_outfile << 4*this->result[i].locus << '\t';
            fs_outfile << this->result[i].T_bf << '\t';
            fs_outfile << this->result[i].T_freq << '\t';
            fs_outfile << this->result[i].T_count << '\t';
            fs_outfile << this->result[i].cvg << '\t';
            fs_outfile << 0 << '\t' << -1 << endl;
        }
        
    }
    
    fs_outfile.close();
}


