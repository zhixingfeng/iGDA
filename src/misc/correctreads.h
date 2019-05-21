//
//  correctreads.h
//  iGDA
//
//  Created by Zhixing Feng on 5/21/19.
//  Copyright © 2019 Zhixing Feng. All rights reserved.
//

#ifndef correctreads_h
#define correctreads_h
#include "io.h"
#include "dist.h"
#include "basic.h"

inline void correct_reads(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range,
                          vector<vector<int> > &encode_data_cd, vector<ReadRange> &reads_range_cd,
                          double min_overlap_rate = 0.75)
{
    if (encode_data.size() != reads_range.size())
        throw runtime_error("correct_reads(): encode_data.size() != reads_range.size()");
    
    size_t genome_size = get_genome_size(reads_range);
    
    vector<bool> temp_array(genome_size*4+3, false);
    
    
    for (auto i = 0; i < encode_data.size(); ++i){
        double max_jaccard = 0;
        int j_max = -1;
        int overlap_start_max = -1;
        int overlap_end_max = -1;
        
        for (auto j = 0; j < encode_data.size(); ++j){
            if (i == j) continue;
            
            // check overlap length
            int overlap_start = reads_range[i].first >= reads_range[j].first? reads_range[i].first : reads_range[j].first;
            int overlap_end = reads_range[i].second <= reads_range[j].second? reads_range[i].second : reads_range[j].second;
            
            int min_overlap = int(min_overlap_rate * (reads_range[i].second - reads_range[i].first + 1));
            if (overlap_end - overlap_start + 1 < min_overlap)
                continue;
            
            // calculate jaccard index
            double cur_jaccard = sim_jaccard(encode_data[i], encode_data[j], reads_range[i], reads_range[j], temp_array, false, min_overlap);
            if (cur_jaccard > max_jaccard){
                max_jaccard = cur_jaccard;
                j_max = j;
                overlap_start_max = overlap_start;
                overlap_end_max = overlap_end;
            }
        }
        
        if (max_jaccard > 0){
            encode_data_cd.push_back(intersect(encode_data[i], encode_data[j_max]));
            reads_range_cd.push_back(pair<int,int>(overlap_start_max, overlap_end_max));
        }else{
            encode_data_cd.push_back(encode_data[i]);
            reads_range_cd.push_back(pair<int,int>(reads_range[i].first, reads_range[i].second));
        }
        
    }
}



#endif /* correctreads_h */
