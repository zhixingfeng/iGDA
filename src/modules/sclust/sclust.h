//
//  sclust.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/7/13.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef sclust_h
#define sclust_h

#include "../../../include/headers.h"
#include "../../misc/misc.h"
#include <thread>
#include <mutex>

class SClust
{
public:
    SClust(){}
    virtual ~SClust(){}
    
    void run(string encode_file, string align_file, string cmpreads_file, 
             string out_file, string tmp_dir, int min_ratio, int min_count, int min_cvg, int n_thread);
    
protected:
    bool run_thread(string cmpreads_file, string out_file, int min_ratio, int min_count, int min_cvg);
    
protected:
    // pileup variants and reads
    vector<vector<int> > pu_var;
    vector<vector<int> > pu_read;
    
    int64_t nreads;
};





#endif /* sclust_h */
