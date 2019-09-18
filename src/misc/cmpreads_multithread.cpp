//
//  cmpreads_multithread.cpp
//  iGDA
//
//  Created by Zhixing Feng on 9/18/19.
//  Copyright © 2019 Zhixing Feng. All rights reserved.
//

#include "cmpreads_multithread.h"

void cmpreads_topn_multithread_core(const vector<vector<int> > &encode_data, const vector<ReadRange> &reads_range, string out_file,
                                           size_t temp_array_size, int start_index, int end_index, int topn, double min_overlap, bool is_rm_single,
                                           bool is_binary, bool is_print_read_id, bool is_condprob, bool is_jaccard)
{
    if (encode_data.size() != reads_range.size())
        throw runtime_error("cmpreads: size of encode_data and reads_range do not match.");
    
    cout << encode_data.size() << endl;
    
    // get the right-most variant location to determing size of template array
    
    vector<int> temp_array(temp_array_size, -1);
    
    // open output file
    FILE *p_out_file = NULL;
    if (is_binary){
        p_out_file = fopen(out_file.c_str(), "wb");
        if (p_out_file==NULL)
            throw runtime_error("fail to open out_file");
    }
    else{
        p_out_file = fopen(out_file.c_str(), "w");
        if (p_out_file==NULL)
            throw runtime_error("fail to open out_file");
        
    }
    
    // pairwise comparison
    for (int i = start_index; i <= end_index; i++){
        if ((i+1)%1000==0)
            cout << i+1 << endl;
        
        // fill the template array by the variants in the ith read
        for (int j = 0; j < encode_data[i].size(); j++)
            temp_array[encode_data[i][j]] = i;
        
        // store matches of the jth read to the ith read
        priority_queue<ReadMatch, vector<ReadMatch>, queue_compare> the_matches;
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (j == i)
                continue;
            
            ReadMatch cur_the_matches;
            
            // get size of overlap of the two reads
            cur_the_matches.start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            cur_the_matches.end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = cur_the_matches.end - cur_the_matches.start + 1;
            
            if (n_overlap < min_overlap * (reads_range[i].second - reads_range[i].first + 1))
                continue;
            
            // get union of the two reads
            int n_union = 0;
            
            for (auto k = 0; k < encode_data[i].size(); ++k){
                if (encode_data[i][k] >= 4*cur_the_matches.start && encode_data[i][k] <= 4*cur_the_matches.end + 3)
                    ++n_union;
            }
            
            for (auto k = 0; k < encode_data[j].size(); ++k){
                if (encode_data[j][k] >= 4*cur_the_matches.start && encode_data[j][k] <= 4*cur_the_matches.end + 3)
                    ++n_union;
            }
            
            // get intersection between two reads
            vector<int> cur_match;
            for (int k = 0; k < encode_data[j].size(); k++){
                if (temp_array[encode_data[j][k]] == i){
                    cur_match.push_back(encode_data[j][k]);
                    --n_union;
                }
            }
            if (n_union < 0)
                throw runtime_error("cmpreads_topn(): n_union < 0");
            
            cur_the_matches.matches = cur_match;
            cur_the_matches.n_overlap = n_overlap;
            
            if (is_jaccard){
                if (n_union > 0)
                    cur_the_matches.match_rate = (double) cur_match.size() / n_union;
                else
                    cur_the_matches.match_rate = 0;
            }else{
                if (is_condprob){
                    if (encode_data[j].size() > 0)
                        cur_the_matches.match_rate = (double) cur_match.size() / encode_data[i].size();
                    else
                        cur_the_matches.match_rate = 0;
                }else{
                    cur_the_matches.match_rate = (double) cur_match.size() / n_overlap;
                }
            }
            
            // keep topn matches
            //the_matches.push(cur_the_matches);
            if (the_matches.size() < topn){
                the_matches.push(cur_the_matches);
            }else{
                if (cur_the_matches.match_rate > the_matches.top().match_rate){
                    the_matches.pop();
                    the_matches.push(cur_the_matches);
                }
            }
        }
        
        // print topn matches
        while(!the_matches.empty()){
            //for (auto i = 0; i < topn; ++i){
            //cout << "the_matches.size() = " << the_matches.size() << endl;
            if (the_matches.empty())
                break;
            ReadMatch tmp_match = the_matches.top();
            // skip matches with size < 2 if is_rm_single is true
            if (is_rm_single){
                if (tmp_match.matches.size() < 2){
                    the_matches.pop();
                    continue;
                }
            }else{
                if (tmp_match.matches.size() == 0){
                    the_matches.pop();
                    continue;
                }
            }
            // print results
            if (is_binary){
                int cur_match_size = (int)tmp_match.matches.size();
                if (is_print_read_id){
                    fwrite(&i, sizeof(int), 1, p_out_file);
                    fwrite(&tmp_match.start, sizeof(int), 1, p_out_file);
                    fwrite(&tmp_match.end, sizeof(int), 1, p_out_file);
                }
                fwrite(&cur_match_size, sizeof(int), 1, p_out_file);
                fwrite(&tmp_match.matches[0], sizeof(int), cur_match_size, p_out_file);
            }else{
                if (is_print_read_id){
                    fprintf(p_out_file, "%d\t", i);
                    fprintf(p_out_file, "%d\t", tmp_match.start);
                    fprintf(p_out_file, "%d\t", tmp_match.end);
                }
                for (int k=0; k<(int)tmp_match.matches.size(); k++)
                    fprintf(p_out_file, "%d,", tmp_match.matches[k]);
                fprintf(p_out_file, "\n");
            }
            the_matches.pop();
        }
    }
    cout << encode_data.size() << endl;
    
    fclose(p_out_file);
}

bool cmpreads_topn_multithread(string encode_file, string align_file, string out_file, int nthread, int topn, double min_overlap,
                                      bool is_rm_single, bool is_binary, bool is_print_read_id, bool is_condprob, bool is_jaccard)
{
    cout << "nthread = " << nthread << endl;
    if (nthread < 1)
        throw runtime_error("cmpreads_topn_multithread(): nthread < 1");
    
    cout << "load encode" << endl;
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    cout << "load reads_range" << endl;
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (encode_data.size() != reads_range.size())
        throw runtime_error("cmpreads: size of encode_data and reads_range do not match.");
    
    cout << encode_data.size() << endl;
    if (encode_data.size() == 0)
        return false;
    
    
    // get the right-most variant location to determing size of template array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    
    // split encode_data and distribute thread
    size_t blocksize_leftover = encode_data.size() % nthread;
    size_t blocksize = encode_data.size() / nthread;
    
    if (blocksize == 0)
        throw runtime_error("cmpreads_topn_multithread(): blocksize == 0");
    
    size_t cur_index = 0;
    vector<pair<size_t, size_t> > block_range;
    while(cur_index < encode_data.size()){
        if (cur_index + blocksize <= encode_data.size()){
            block_range.push_back(pair<size_t, size_t>(cur_index, cur_index + blocksize - 1));
        }else{
            if (blocksize_leftover == 0)
                throw runtime_error("cmpreads_topn_multithread(): blocksize_leftover == 0");
            block_range.back().second += blocksize_leftover;
        }
        cur_index += blocksize;
    }
    
    // multithread run
    
    vector<thread> threads;
    for (auto i = 0; i < block_range.size(); ++i){
        //cout << block_range[i].first << "\t" << block_range[i].second << endl;
        if (block_range[i].first < 0 || block_range[i].second >= encode_data.size())
                throw runtime_error("cmpreads_topn_multithread(): block_range[i].first < 0 || block_range[i].second >= encode_data.size()");
        //string cur_out_file = out_file + ".part_" + to_string(block_range[i].first) + "_" + to_string(block_range[i].second);
        string cur_out_file = out_file + ".part_" + to_string(i);
        threads.push_back(thread(cmpreads_topn_multithread_core, std::ref(encode_data), std::ref(reads_range), cur_out_file, temp_array_size,
                                 (int)block_range[i].first, (int)block_range[i].second, topn, min_overlap, is_rm_single, is_binary,
                                 is_print_read_id, is_condprob, is_jaccard));
    }
     
    for (auto i = 0; i < threads.size(); ++i)
        threads[i].join();
    
    return true;
  
}


