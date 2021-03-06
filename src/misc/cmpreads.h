//
//  cmpreads.h
//  iGDA
//
//  Created by Zhixing Feng on 16/9/7.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.

//  pairwise compare encoded reads to find consistently occuring variants

#ifndef cmpreads_h
#define cmpreads_h

#include "io.h"


/*inline bool operator==(const ReadMatch& lhs, const ReadMatch& rhs){ return lhs.match_rate == rhs.match_rate; }
inline bool operator<=(const ReadMatch& lhs, const ReadMatch& rhs){ return lhs.match_rate <= rhs.match_rate; }
inline bool operator<(const ReadMatch& lhs, const ReadMatch& rhs){ return lhs.match_rate < rhs.match_rate; }
inline bool operator>=(const ReadMatch& lhs, const ReadMatch& rhs){ return lhs.match_rate >= rhs.match_rate; }
inline bool operator>(const ReadMatch& lhs, const ReadMatch& rhs){ return lhs.match_rate > rhs.match_rate; }
*/
struct queue_compare
{
    bool operator()(const ReadMatch& l, const ReadMatch& r)
    {
        if(l.match_rate != r.match_rate)
            return l.match_rate > r.match_rate;
        
        if (l.matches.size() != r.matches.size())
            return l.matches.size() > r.matches.size();
        
        if (l.matches.size() != 0)
            return l.matches[0] > r.matches[0];
        
        return true;
    }
};


inline bool readmatch_compare(const ReadMatch& l, const ReadMatch& r)
{
    return l.end > r.end;
}

inline bool group_cmpreadsdiff(const vector<ReadMatch> &cmpreadsdiff_data, vector<vector<ReadMatch> > &cmpreadsdiff_data_grouped, bool is_sort = true)
{
    // get maximal code
    int max_code = 0;
    for (size_t i = 0; i < cmpreadsdiff_data.size(); ++i)
        for (size_t j = 0; j < cmpreadsdiff_data[i].matches.size(); ++j)
            if (cmpreadsdiff_data[i].matches[j] > max_code)
                max_code = cmpreadsdiff_data[i].matches[j];
        
    // group cmpreadsdiff
    cmpreadsdiff_data_grouped.clear();
    cmpreadsdiff_data_grouped.resize(max_code + 1);
    
    for (size_t i = 0; i < cmpreadsdiff_data.size(); ++i){
        if (cmpreadsdiff_data[i].matches.size() > 0){
            cmpreadsdiff_data_grouped[cmpreadsdiff_data[i].matches[0]].push_back(cmpreadsdiff_data[i]);
        }
    }
    
    // sort cmpreadsdiff_data_grouped
    if (is_sort){
        for (size_t i = 0; i < cmpreadsdiff_data_grouped.size(); ++i){
            if (cmpreadsdiff_data_grouped[i].size() > 0){
                vector<ReadMatch> buf = cmpreadsdiff_data_grouped[i];
                sort(buf.begin(), buf.end(), readmatch_compare);
                cmpreadsdiff_data_grouped[i] = buf;
            }
        }
    }
    
    return true;
}

// compare reads and use top n as candidates (include difference between reads, read IDs will be added)
inline bool cmpreads_topn_diff(string encode_file, string align_file, string out_file, int topn = 20, double min_overlap = 0.5)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (encode_data.size() != reads_range.size())
        throw runtime_error("cmpreads: size of encode_data and reads_range do not match.");
    
    cout << encode_data.size() << endl;
    
    // get the right-most variant location to determing size of template array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // open output file
    FILE *p_out_file = fopen(out_file.c_str(), "wb");
    if (p_out_file==NULL)
        throw runtime_error("fail to open out_file");
    
    // pairwise comparison
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0) cout << i+1 << endl;
        
        // fill the template array by the variants in the ith read
        for (int j = 0; j < encode_data[i].size(); j++)
            temp_array[encode_data[i][j]] = i;
        
        // store matches of the jth read to the ith read (only store topn using priority_queue)
        priority_queue<ReadMatch, vector<ReadMatch>, queue_compare> the_matches;
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (j == i)
                continue;
            
            // get size of overlap of the two reads
            int overlap_start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            int overlap_end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = overlap_end - overlap_start + 1;
            
            if (n_overlap < min_overlap * (reads_range[i].second - reads_range[i].first + 1))
                continue;
            
            
            // get intersection and difference between two reads
            vector<int> cur_match;
            vector<int> cur_diff;
            for (int k = 0; k < encode_data[j].size(); k++){
                if (temp_array[encode_data[j][k]] == i){
                    cur_match.push_back(encode_data[j][k]);
                }else{
                    if (encode_data[j][k] >= 4*reads_range[i].first && encode_data[j][k] <= 4*reads_range[i].second + 3)
                        cur_diff.push_back(encode_data[j][k]);
                }
            }
            
            // calculate match rate
            double cur_match_rate = 0;
            if (encode_data[j].size() > 0)
                cur_match_rate = (double) cur_match.size() / encode_data[i].size();
            
            // keep topn matches
            if (the_matches.size() < topn){
                the_matches.push(ReadMatch(cur_diff, cur_match, cur_match_rate, n_overlap, i, overlap_start, overlap_end));
            }else{
                if (cur_match_rate > the_matches.top().match_rate){
                    the_matches.pop();
                    the_matches.push(ReadMatch(cur_diff, cur_match, cur_match_rate, n_overlap, i, overlap_start, overlap_end));
                }
            }
            
        }
        
        // print topn matches
        while(!the_matches.empty()){
            ReadMatch tmp_match = the_matches.top();
            if (tmp_match.matches.size() >= 2){
                int cur_match_size = (int)tmp_match.matches.size();
                int cur_diff_size = (int)tmp_match.diff.size();
                
                fwrite(&i, sizeof(int), 1, p_out_file);
                fwrite(&tmp_match.start, sizeof(int), 1, p_out_file);
                fwrite(&tmp_match.end, sizeof(int), 1, p_out_file);
                
                fwrite(&cur_match_size, sizeof(int), 1, p_out_file);
                fwrite(&tmp_match.matches[0], sizeof(int), cur_match_size, p_out_file);
                
                fwrite(&cur_diff_size, sizeof(int), 1, p_out_file);
                fwrite(&tmp_match.diff[0], sizeof(int), cur_diff_size, p_out_file);
            }
            the_matches.pop();
        }
    }
    
    fclose(p_out_file);
    return true;
}

// compare reads and use top n as candidates (read data from files)
inline bool cmpreads_topn(string encode_file, string align_file, string out_file, int topn = 20, double min_overlap = 0.5,
                          bool is_rm_single=true, bool is_binary=true, bool is_print_read_id=false, bool is_condprob=true, bool is_jaccard = true)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (encode_data.size() != reads_range.size())
        throw runtime_error("cmpreads: size of encode_data and reads_range do not match.");
    
    cout << encode_data.size() << endl;
    
    // get the right-most variant location to determing size of template array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
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
    for (int i=0; i<(int)encode_data.size(); i++){
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
    
    return true;
}


// compare reads and use top n as candidates (read data from files) (legacy version)
inline bool cmpreads_topn_legacy(string encode_file, string align_file, string out_file, int topn = 20, double min_overlap = 0.5,
                     bool is_rm_single=true, bool is_binary=true, bool is_print_read_id=false, bool is_condprob=true)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (encode_data.size() != reads_range.size())
        throw runtime_error("cmpreads: size of encode_data and reads_range do not match.");
    
    cout << encode_data.size() << endl;
    
    // get the right-most variant location to determing size of template array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
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
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0) cout << i+1 << endl;
        
        // fill the template array by the variants in the ith read
        for (int j = 0; j < encode_data[i].size(); j++)
            temp_array[encode_data[i][j]] = i;
        
        // store matches of the jth read to the ith read
        vector<ReadMatch> the_matches (encode_data.size(), ReadMatch());
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (j == i)
                continue;
            
            // get size of overlap of the two reads
            the_matches[j].start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            the_matches[j].end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = the_matches[j].end - the_matches[j].start + 1;
           
            if (n_overlap < min_overlap * (reads_range[i].second - reads_range[i].first + 1))
                continue;
            
            
            // get intersection between two reads
            vector<int> cur_match;
            for (int k = 0; k < encode_data[j].size(); k++)
                if (temp_array[encode_data[j][k]] == i)
                    cur_match.push_back(encode_data[j][k]);
            
            
            the_matches[j].matches = cur_match;
            the_matches[j].n_overlap = n_overlap;
            
            if (is_condprob){
                if (encode_data[j].size() > 0)
                    the_matches[j].match_rate = (double) cur_match.size() / encode_data[i].size();
                else    
                    the_matches[j].match_rate = 0;
            }else{
                the_matches[j].match_rate = (double) cur_match.size() / n_overlap;
            }
        }
        
        // sort the_matches according to match_rate
        stable_sort(the_matches.begin(), the_matches.end(), [](const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate > dr.match_rate;});
        
        // print topn matches
        int cur_size = topn <= the_matches.size() ? topn : (int)the_matches.size();
        for (int j=0; j<cur_size; j++){
            // skip matches with size < 2 if is_rm_single is true
            if (is_rm_single){
                if (the_matches[j].matches.size() < 2)
                    continue;
            }else{
                if (the_matches[j].matches.size() == 0)
                    continue;
            }
            // print results
            if (is_binary){
                int cur_match_size = (int)the_matches[j].matches.size();
                if (is_print_read_id){
                    fwrite(&i, sizeof(int), 1, p_out_file);
                    fwrite(&the_matches[j].start, sizeof(int), 1, p_out_file);
                    fwrite(&the_matches[j].end, sizeof(int), 1, p_out_file);
                }
                fwrite(&cur_match_size, sizeof(int), 1, p_out_file);
                fwrite(&the_matches[j].matches[0], sizeof(int), cur_match_size, p_out_file);
            }else{
                if (is_print_read_id){
                    fprintf(p_out_file, "%d\t", i);
                    fprintf(p_out_file, "%d\t", the_matches[j].start);
                    fprintf(p_out_file, "%d\t", the_matches[j].end);
                }
                for (int k=0; k<(int)the_matches[j].matches.size(); k++)
                    fprintf(p_out_file, "%d,", the_matches[j].matches[k]);
                fprintf(p_out_file, "\n");
            }

            
        }
        
    }
    cout << encode_data.size() << endl;
    
    fclose(p_out_file);
    
    return true;
}



// convert binary cmpreadsfile to text file
inline void cmpreads_bin2txt(string cmpreads_binfile, string cmpreads_txtfile, bool is_read_id = false, bool is_diff = false)
{
    // open binary file
    FILE *p_binfile = fopen(cmpreads_binfile.c_str(), "rb");
    if (p_binfile == NULL)
        runtime_error("fail to open cmpreads_binfile");
    
    // open text file to be written
    FILE *p_txtfile = fopen(cmpreads_txtfile.c_str(), "wb");
    if (p_txtfile == NULL)
        runtime_error("fail to open cmpreads_txtfile");   
    
    // scan binary file and convert
    if (is_diff){
        while(1){
            int cand_loci_size;
            int cand_loci_diff_size;
            int read_id;
            int start;
            int end;
            fread(&read_id, sizeof(int), 1, p_binfile);
            fread(&start, sizeof(int), 1, p_binfile);
            fread(&end, sizeof(int), 1, p_binfile);
            
            fread(&cand_loci_size, sizeof(int), 1, p_binfile);
            vector<int> cand_loci(cand_loci_size,-1);
            fread(&cand_loci[0], sizeof(int), cand_loci_size, p_binfile);
            
            fread(&cand_loci_diff_size, sizeof(int), 1, p_binfile);
            vector<int> cand_loci_diff(cand_loci_diff_size,-1);
            fread(&cand_loci_diff[0], sizeof(int), cand_loci_diff_size, p_binfile);
            
            
            if (feof(p_binfile))
                break;
            
            fprintf(p_txtfile, "%d\t", read_id);
            fprintf(p_txtfile, "%d\t", start);
            fprintf(p_txtfile, "%d\t", end);
            
            for (int i = 0; i < cand_loci_size; i++)
                fprintf(p_txtfile, "%d,", cand_loci[i]);
            
            fprintf(p_txtfile, "\t");
            
            for (int i = 0; i < cand_loci_diff_size; i++)
                fprintf(p_txtfile, "%d,", cand_loci_diff[i]);
            fprintf(p_txtfile, "\n");
        }
    }else{
        if (is_read_id){
            while(1){
                int cand_loci_size;
                int read_id;
                int start;
                int end;
                fread(&read_id, sizeof(int), 1, p_binfile);
                fread(&start, sizeof(int), 1, p_binfile);
                fread(&end, sizeof(int), 1, p_binfile);
                fread(&cand_loci_size, sizeof(int), 1, p_binfile);
                vector<int> cand_loci(cand_loci_size,-1);
                fread(&cand_loci[0], sizeof(int), cand_loci_size, p_binfile);
                if (feof(p_binfile))
                    break;
                
                fprintf(p_txtfile, "%d\t", read_id);
                fprintf(p_txtfile, "%d\t", start);
                fprintf(p_txtfile, "%d\t", end);
                for (int i = 0; i < cand_loci_size; i++)
                    fprintf(p_txtfile, "%d,", cand_loci[i]);
                fprintf(p_txtfile, "\n");
            }
        }else{
            while(1){
                int cand_loci_size;
                fread(&cand_loci_size, sizeof(int), 1, p_binfile);
                vector<int> cand_loci(cand_loci_size,-1);
                fread(&cand_loci[0], sizeof(int), cand_loci_size, p_binfile);
                if (feof(p_binfile))
                    break;
                for (int i = 0; i < cand_loci_size; i++)
                    fprintf(p_txtfile, "%d,", cand_loci[i]);
                fprintf(p_txtfile, "\n");
            }   
        }
    }
    fclose(p_binfile);
    fclose(p_txtfile);
    
}

// convert text cmpreadsfile to binary file  
inline void cmpreads_txt2bin(string cmpreads_txtfile, string cmpreads_binfile, bool is_read_id=false)
{
    // open text file
    ifstream p_txtfile;
    open_infile(p_txtfile, cmpreads_txtfile);
    // open binary file
    FILE *p_binfile = fopen(cmpreads_binfile.c_str(), "wb");
    if (p_binfile == NULL)
        runtime_error("fail to open cmpreads_binfile");
    
    while(1){
        // read text file
        string buf;
        getline(p_txtfile, buf);
        if (p_txtfile.eof())
            break;
        
        if (is_read_id){
            vector<string> cur_data = split(buf, '\t');
            if(cur_data.size() != 4)
                throw runtime_error("incorrect format in cmpreads_txtfile");
            
            int read_id = stoi(cur_data[0]);
            int start = stoi(cur_data[1]);
            int end = stoi(cur_data[2]);
            vector<int> cand_loci = split_int(cur_data[3], ',');
            int cand_loci_size = (int) cand_loci.size();
        
            // write binary file
            fwrite(&read_id, sizeof(int), 1, p_binfile);
            fwrite(&start, sizeof(int), 1, p_binfile);
            fwrite(&end, sizeof(int), 1, p_binfile);
            fwrite(&cand_loci_size, sizeof(int), 1, p_binfile);
            fwrite(&cand_loci[0], sizeof(int), cand_loci_size, p_binfile);
        }else{
            vector<int> cand_loci = split_int(buf, ',');
            int cand_loci_size = (int) cand_loci.size();
            
            // write binary file
            fwrite(&cand_loci_size, sizeof(int), 1, p_binfile);
            fwrite(&cand_loci[0], sizeof(int), cand_loci_size, p_binfile);
        }
    }
    p_txtfile.close();
    fclose(p_binfile);
}

inline map<string, vector<int64_t> > cmpreads_split(string cmpreads_binfile, string prefix, int n_parts)
{
    map<string, vector<int64_t> > rl;
    // open binary file
    FILE *p_binfile = fopen(cmpreads_binfile.c_str(), "rb");
    if (p_binfile == NULL)
        runtime_error("fail to open cmpreads_binfile");
    
    // get number of candidates
    int64_t n_cand = 0;
    while(1){
        int cand_loci_size;
        fread(&cand_loci_size, sizeof(int), 1, p_binfile);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_binfile);
        if (feof(p_binfile))
            break;
        n_cand++;
    }
    fclose(p_binfile);
    
    cout << "n_cand: " << n_cand << endl;
    
    
    // calculate size of each part
    vector<int64_t> part_sizes(n_parts, (int64_t) n_cand / n_parts);
    part_sizes[0] = (int64_t) n_cand / n_parts + (int64_t) n_cand % n_parts;
    
    // get start location of the file for each part
    p_binfile = fopen(cmpreads_binfile.c_str(), "rb");
    if (p_binfile == NULL)
        runtime_error("fail to open cmpreads_binfile");
    
    for (int i=0; i<(int)part_sizes.size(); i++){
        string outfile = prefix + "_" + to_string(i);
        FILE *p_outfile = fopen(outfile.c_str(), "wb");
        if (p_outfile == NULL)
            runtime_error("fail to open outfile");
        for (int j=0; j<part_sizes[i]; j++){
            // read a candidate
            int cand_loci_size;
            fread(&cand_loci_size, sizeof(int), 1, p_binfile);
            vector<int> cand_loci(cand_loci_size,-1);
            fread(&cand_loci[0], sizeof(int), cand_loci_size, p_binfile);
            if (feof(p_binfile))
                break;
            
            // write the candidate
            fwrite(&cand_loci_size, sizeof(int), 1, p_outfile);
            fwrite(&cand_loci[0], sizeof(int), cand_loci_size, p_outfile);
        }
        
        fclose(p_outfile);
    }
    
    fclose(p_binfile);
    return rl;
}
#endif /* cmpreads_h */








