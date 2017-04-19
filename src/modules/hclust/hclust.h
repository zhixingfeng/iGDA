//
//  hclust.hpp
//  iGDA
//
//  Created by Zhixing Feng on 17/4/18.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#ifndef hclust_h
#define hclust_h

#include "../../../include/headers.h"
#include "../../misc/io.h"
#include "../aligncoder/aligncodersnv.h"

class HClust{
    
public:
    HClust():ptr_aligncoder(NULL) {}
    HClust(AlignCoder *a_ptr_aligncoder) {ptr_aligncoder = a_ptr_aligncoder;}
    virtual ~HClust() {}
    void setAlignCoder(AlignCoder *a_ptr_aligncoder) {ptr_aligncoder = a_ptr_aligncoder;}
    
    // mask variants outside of region_file, which contains a list of loci (by default 1-based)
    void mask(string encode_file, string region_file, string out_file, bool is_0_based = false);
    
    // caculate distance
    void dist(string encode_file, string align_file, string out_file, bool is_nmiss=false);
    
protected:
    AlignCoder *ptr_aligncoder;
    
    
};


#endif /* hclust_h */
