//
//  cmpreads_multithread.h
//  iGDA
//
//  Created by Zhixing Feng on 9/16/19.
//  Copyright © 2019 Zhixing Feng. All rights reserved.
//

#ifndef cmpreads_multithread_h
#define cmpreads_multithread_h
#include "cmpreads.h"
#include <thread>
void cmpreads_topn_multithread_core(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range, string out_file,
                                    size_t temp_array_size, int start_index, int end_index, int topn = 20, double min_overlap = 0.5, bool is_rm_single=true,
                                    bool is_binary=true, bool is_print_read_id=false, bool is_condprob=true, bool is_jaccard = true);

bool cmpreads_topn_multithread(string encode_file, string align_file, string out_file, int nthread = 1, int topn = 20, double min_overlap = 0.5,
                               bool is_rm_single=true, bool is_binary=true, bool is_print_read_id=false, bool is_condprob=true, bool is_jaccard = true);

#endif /* cmpreads_multithread_h */
