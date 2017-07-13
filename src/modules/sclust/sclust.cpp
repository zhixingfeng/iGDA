//
//  sclust.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/7/13.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "sclust.h"

void SClust::run(string encode_file, string align_file, string cmpreads_file, 
         string out_file, string tmp_dir, int min_cvg, int n_thread)
{
    // pileup reads 
    int64_t nreads_pu_var, nreads_pu_read;
    vector<vector<int> > pu_var = pileup_var(encode_file, nreads_pu_var);
    vector<vector<int> > pu_read = pileup_reads(align_file, nreads_pu_read);
    
    // run the subpace clustering algorithm
    if (n_thread == 1){
        // single thread
        
    }else{
        // multiple threads
    }
}
