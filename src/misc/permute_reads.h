//
//  permute_reads.hpp
//  iGDA
//
//  Created by Zhixing Feng on 1/28/20.
//  Copyright © 2020 Zhixing Feng. All rights reserved.
//

#ifndef permute_reads_hpp
#define permute_reads_hpp

#include "../../include/headers.h"
#include "./statistics.h"

void permute_encodefile(string m5_file, string pileup_file, string outfile, int seed = 18473);
void get_condprob_threshold(string dforest_permuted_file, string pileup_file);


#endif /* permute_reads_hpp */
