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

struct ReadMatch
{
    ReadMatch():match_rate(0),n_overlap(0){}
    vector<int> matches;
    double match_rate;
    int n_overlap;
    int start;
    int end;
};

/*inline bool operator < (const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate < dr.match_rate;}
inline bool operator <= (const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate <= dr.match_rate;}
inline bool operator > (const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate > dr.match_rate;}
inline bool operator >= (const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate >= dr.match_rate;}
inline bool operator == (const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate == dr.match_rate;}
*/

inline bool cmpreads_topn(string encode_file, string align_file, string out_file, int topn = 10, double min_overlap = 0.25,
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
            //if (n_overlap < min_overlap * (reads_range[i].second - reads_range[i].first + 1) && 
            //    n_overlap < min_overlap * (reads_range[j].second - reads_range[j].first + 1))
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
                    the_matches[j].match_rate = (double) cur_match.size() / encode_data[j].size();
                else    
                    the_matches[j].match_rate = 0;
            }else{
                the_matches[j].match_rate = (double) cur_match.size() / n_overlap;
            }
        }
        
        // sort the_matches according to match_rate
        sort(the_matches.begin(), the_matches.end(), [](const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate > dr.match_rate;});
        
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



/*inline bool cmpreads(string encode_file, string align_file, string out_file, double min_overlap = 0.25,
                     bool is_rm_single=true, bool is_binary=true)
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
    
    // pairwise comparison
    FILE *p_out_binfile = NULL;
    FILE *p_out_textfile = NULL;
    if (is_binary){
        p_out_binfile = fopen(out_file.c_str(), "wb");
        if (p_out_binfile==NULL)
            throw runtime_error("fail to open out_file");
    }
    else{
        p_out_textfile = fopen(out_file.c_str(), "w");
        if (p_out_textfile==NULL)
            throw runtime_error("fail to open out_file");

    }
    
    
    for (int i=0; i<(int)(encode_data.size()-1); i++){
        if ((i+1)%1000==0) cout << i+1 << endl;
        
        // fill the template array by the variants in the ith read
        for (int j = 0; j < encode_data[i].size(); j++)
            temp_array[encode_data[i][j]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=i+1; j<(int)encode_data.size(); j++){
            //if (reads_range[i].first > reads_range[j].second || reads_range[j].first > reads_range[i].second)
            //    continue;
            
            int n_overlap = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second -
                            reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first + 1;
            if (n_overlap < min_overlap * (reads_range[i].second - reads_range[i].first + 1) && n_overlap < min_overlap * (reads_range[j].second - reads_range[j].first + 1))
                continue;
            vector<int> cur_match;
            for (int k = 0; k < encode_data[j].size(); k++)
                if (temp_array[encode_data[j][k]] == i)
                    cur_match.push_back(encode_data[j][k]);
            
            if (is_rm_single){
                if (cur_match.size() < 2)
                    continue;
            }else{
                if (cur_match.size() == 0)
                    continue;
            }
            if (is_binary){
                int cur_match_size = (int)cur_match.size();
                fwrite(&cur_match_size, sizeof(int), 1, p_out_binfile);
                fwrite(&cur_match[0], sizeof(int), cur_match_size, p_out_binfile);
            }else{
                for (int k=0; k<(int)cur_match.size(); k++)
                    fprintf(p_out_textfile, "%d,", cur_match[k]);
                fprintf(p_out_textfile, "\n");
            }
        }
                
    }
    cout << encode_data.size() << endl;
    
    if (is_binary)
        fclose(p_out_binfile);
    else
        fclose(p_out_textfile);
    
    return true;
}*/

// convert binary cmpreadsfile to text file
inline void cmpreads_bin2txt(string cmpreads_binfile, string cmpreads_txtfile, bool is_read_id=false)
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








