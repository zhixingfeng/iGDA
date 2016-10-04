//
//  pileup.h
//  iGDA
//
//  Created by Zhixing Feng on 16/10/3.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef pileup_h
#define pileup_h

#include "./cmpreads.h"

// coverage should NOT exceed range of int

inline vector<int> pileup(string encode_file)
{
    // load encode data
    cout << "load encode" << endl;
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // calculate size of pileup vector
    cout << "calculate size of pileup vector" << endl;
    int pu_size = 0;
    for (int i=0; i<(int)encode_data.size(); i++)
        for (int j=0; j<(int)encode_data[i].size(); j++)
            pu_size = encode_data[i][j] > pu_size ? encode_data[i][j] : pu_size;
    
    // pileup
    cout << "pileup" << endl;
    vector<int> pu(pu_size, 0);
    for (int i=0; i<(int)encode_data.size(); i++)
        for (int j=0; j<(int)encode_data[i].size(); j++)
            pu[encode_data[i][j]-1]++;
    
    return pu;
}


#endif /* pileup_h */
