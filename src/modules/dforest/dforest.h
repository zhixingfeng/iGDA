//
//  dforest.h
//  iGDA
//
//  Created by Zhixing Feng on 16/12/1.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//
//  To implement data-based forest method for bump hunting, we demond large amount of computation as we have to train p*n^2 models,
//  where p is average shared variants of two reads and n is number of reads. We need fast yet low memory cost ways to train models.
//  For large datasets, read ID pileup file might not be able to fit into memory, so we use two different way to implement: first is
//  the naive way, and the second is prefix compression (memory efficient but might be slower).

//  Input is pileup data, i.e Pileup variant and pileup read ID.
//  Parameters: minimal reads supporting a partition, maximal depth of a tree.


#ifndef dforest_h
#define dforest_h

#include "../../../include/headers.h"

// define final result.
struct Result{
    double bf;
    double prop;
    vector<int> candidate_id;
    vector<int> link_loci;
};







#endif /* dforest_h */
