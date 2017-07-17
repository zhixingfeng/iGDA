//
//  sclust.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/7/13.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "sclust.h"
mutex mtx_sclust;
void SClust::run(string encode_file, string align_file, string cmpreads_file, 
                 string out_file, string tmp_dir, int min_ratio, int min_count, int min_cvg, int n_thread)
{
    // pileup reads 
    int64_t nreads_pu_var, nreads_pu_read;
    pu_var = pileup_var(encode_file, nreads_pu_var);
    pu_read = pileup_reads(align_file, nreads_pu_read);
    
    // run the subpace clustering algorithm
    if (n_thread == 1){
        // single thread
        run_thread(cmpreads_file, out_file, min_ratio, min_cvg, min_cvg);
    }else{
        // multiple threads
        
    }
}

bool SClust::run_thread(string cmpreads_file, string out_file, int min_ratio, int min_count, int min_cvg){
    
    return true;
}
