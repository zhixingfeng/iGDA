//
//  dforest.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/12/1.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "dforest.h"

void DForest::filter(string dforest_file, string out_file, double minfreq)
{
    ifstream fs_infile;
    ofstream fs_outfile;
    open_infile(fs_infile, dforest_file);
    open_outfile(fs_outfile, out_file);
    
    int64_t k = 1;
    while(!fs_infile.eof()){
        string buf;
        getline(fs_infile, buf);
        if(fs_infile.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if ((int) buf_vec.size() != 7)
            throw runtime_error("incorrect format in " + dforest_file);
        
        if (k%100000 == 0)
            cout << k << endl;
        if (stod(buf_vec[2]) >= minfreq)
            fs_outfile << buf << endl;
        k++;
    }
    
    fs_infile.close();
    fs_outfile.close();
}
