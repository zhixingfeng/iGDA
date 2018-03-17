//
//  assemble.cpp
//  iGDA
//
//  Created by Zhixing Feng on 17/8/29.
//  Copyright © 2017年 Zhixing Feng. All rights reserved.
//

#include "assembler.h"


void Assembler::get_variants(string dforest_file, string out_file, double min_condprob)
{
    AlignCoderSNV aligncodersnv;
    ifstream fs_dforeset_file;
    open_infile(fs_dforeset_file, dforest_file);
    ofstream fs_out_file;
    open_outfile(fs_out_file, out_file);
    while(1){
        string buf;
        getline(fs_dforeset_file, buf);
        if (fs_dforeset_file.eof())
            break;
        vector<string> buf_vec = split(buf,'\t');
        if (buf_vec.size()!=7)
            throw runtime_error("incorrect format in get_variants()");
        
        int code = stoi(buf_vec[0]);
        double condprob = stod(buf_vec[2]);
        if (condprob >= min_condprob){
            pair<int, char> rl_decode = aligncodersnv.decode(code);
            int locus = rl_decode.first;
            char base = rl_decode.second;
            fs_out_file << locus << '\t' << base << '\t' << buf << endl;
        }
        
    }
    fs_out_file.close();
    fs_dforeset_file.close();
    fs_out_file.close();
}

void Assembler::reduce_dim(string encode_file, string var_file, string out_file)
{
    // scan var_file to get maximal code
    ifstream fs_var_file;
    open_infile(fs_var_file, var_file);
    int max_code = 0;
    while (1){
        string buf;
        getline(fs_var_file, buf);
        if (fs_var_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in var_file");
        int cur_code = stoi(buf_vec[2]);
        if ( cur_code > max_code)
            max_code = cur_code;
    }
    fs_var_file.close();

    // fill template by scaning var_file for the second time
    vector<bool> temp_code(max_code + 1, false);
    open_infile(fs_var_file, var_file);
    while (1){
        string buf;
        getline(fs_var_file, buf);
        if (fs_var_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in var_file");
        int cur_code = stoi(buf_vec[2]);
        temp_code[cur_code] = true;
    }
    fs_var_file.close();
    
    // scan encode_file
    ifstream fs_encode_file;
    ofstream fs_out_file;
    open_infile(fs_encode_file, encode_file);
    open_outfile(fs_out_file, out_file);
    while (1){
        string buf;
        getline(fs_encode_file, buf);
        if (fs_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (int i=0; i<(int)buf_vec.size(); ++i){
            if (buf_vec[i] <= max_code){
                if (temp_code[buf_vec[i]]){
                    fs_out_file << buf_vec[i] << '\t';
                }
            }
        }
        fs_out_file << endl;
    }
    fs_var_file.close();
    fs_out_file.close();


}

void Assembler::dist(string encode_file, string align_file, string out_file)
{
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_file (m5 format)
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file, 'm');
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data and reads_range are different");
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0)
            cout << "processed " << i+1 << " reads" <<endl;
        
        // fill the template array by the variants in the ith read
        for (int k = 0; k < encode_data[i].size(); k++)
            temp_array[encode_data[i][k]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (i==j)
                continue;
            int overlap_start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            int overlap_end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = overlap_end - overlap_start + 1;
            //int n_overlap = (reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second) -
            //     (reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first) + 1;
            
            if (n_overlap <= 0)
                continue;
            
            // scan read i, only retain variants in [overlap_start, overlap_end]
            int n_miss = 0;
            for (int k = 0; k < encode_data[i].size(); k++){
                if (encode_data[i][k] >= overlap_start*4 &&
                    encode_data[i][k] <= overlap_end*4 + 3){
                    ++n_miss;
                }
            }
            
            // scan read j, only retain variants in [overlap_start, overlap_end]
            for (int k = 0; k < encode_data[j].size(); k++){
                if (encode_data[j][k] >= overlap_start*4 &&
                    encode_data[j][k] <= overlap_end*4 + 3){
                    if (temp_array[encode_data[j][k]] == i)
                        --n_miss;
                    else
                        ++n_miss;
                    
                }
            }
            
            // print results
            fs_out << i <<',' << j << ',' <<(double)n_miss / n_overlap << ',' << n_miss << ',' << n_overlap << endl;
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();
}


void Assembler::dist_rdim(string encode_file, string align_file, string var_file, string out_file)
{
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_file (m5 format)
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file, 'm');
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data and reads_range are different");
    
    // get maximal genome posotion
    int genome_size = 0;
    for (int i=0; i<(int)reads_range.size(); ++i)
        genome_size = genome_size < reads_range[i].second ? reads_range[i].second : genome_size;
    
    ++genome_size;
    
    // cumulated cdf of var
    ifstream fs_var_file;
    open_infile(fs_var_file, var_file);
    vector<int> cdf_var(genome_size, 0);
    int cur_count = 0;
    int pre_pos = 0;
    while (1){
        string buf;
        getline(fs_var_file, buf);
        if (fs_var_file.eof())
            break;
        vector<string> buf_vec = split(buf, '\t');
        if (buf_vec.size()!=9)
            throw runtime_error("incorrect format in var_file");
        int cur_pos = stod(buf_vec[0]);
        ++cur_count;
        for (int i=pre_pos; i<=cur_pos; ++i)
            cdf_var[i] = cur_count;
        pre_pos = cur_pos + 1;
    }
    fs_var_file.close();
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0)
            cout << "processed " << i+1 << " reads" <<endl;
        
        // fill the template array by the variants in the ith read
        for (int k = 0; k < encode_data[i].size(); k++)
            temp_array[encode_data[i][k]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (i==j)
                continue;
            int overlap_start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            int overlap_end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = overlap_end - overlap_start + 1;
           
            // calculate number of variants in [overlap_start, overlap_end]
            int n_var = cdf_var[overlap_end] - cdf_var[overlap_start];

            if (n_overlap <= 0 || cdf_var[overlap_end] - cdf_var[overlap_start] <=0)
                continue;
            
            // scan read i, only retain variants in [overlap_start, overlap_end]
            int n_miss = 0;
            for (int k = 0; k < encode_data[i].size(); k++){
                if (encode_data[i][k] >= overlap_start*4 &&
                    encode_data[i][k] <= overlap_end*4 + 3){
                    ++n_miss;
                }
            }
            
            // scan read j, only retain variants in [overlap_start, overlap_end]
            for (int k = 0; k < encode_data[j].size(); k++){
                if (encode_data[j][k] >= overlap_start*4 &&
                    encode_data[j][k] <= overlap_end*4 + 3){
                    if (temp_array[encode_data[j][k]] == i)
                        --n_miss;
                    else
                        ++n_miss;
                    
                }
            }
            
            // print results
            fs_out << i <<',' << j << ',' <<(double)n_miss / n_var << ',' << n_miss << ',' << n_overlap <<',' << n_var<< endl;
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();

}



void Assembler::jaccard_index(string encode_file, string align_file, string out_file)
{
    // load encode_file
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load align_file (m5 format)
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file, 'm');
    if (encode_data.size() != reads_range.size())
        throw runtime_error("size of encode_data and reads_range are different");
    
    // scan encode_data to determine size of temp_array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size()-1; i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // pairwise comparison
    ofstream fs_out;
    open_outfile(fs_out, out_file);
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0)
            cout << "processed " << i+1 << " reads" <<endl;
        
        // fill the template array by the variants in the ith read
        for (int k = 0; k < encode_data[i].size(); k++)
            temp_array[encode_data[i][k]] = i;
        
        
        // compare other reads to cur_variant
        for (int j=i+1; j<(int)encode_data.size(); j++){
            //if (i==j) continue;
            
            int start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            int end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = end - start + 1;
            
            int code_start = 4*start;
            int code_end = 4*end + 3;
            
            if (n_overlap <= 0)
                continue;
            
            int n_intersect = 0;
            int n_union = 0;
            for (int k = 0; k < encode_data[i].size(); k++){
                if (encode_data[i][k]>=code_start && encode_data[i][k]<=code_end)
                    ++n_union;
            }
            
            for (int k = 0; k < encode_data[j].size(); k++){
                if (encode_data[j][k] < code_start || encode_data[j][k] > code_end)
                    continue;
                if (temp_array[encode_data[j][k]] == i)
                    ++n_intersect;
                else
                    ++n_union;
            }
            
            double jaccard_index = n_union==0 ? 0 : (double)n_intersect / n_union;
            
            fs_out << i <<',' << j << ',' << jaccard_index << ',';
            fs_out << n_intersect << ',' << n_union << ',' << n_overlap << ',';
            fs_out << start << ',' << end << ',' << code_start << ',' << code_end << endl;
            
            fs_out << j <<',' << i << ',' << jaccard_index << ',';
            fs_out << n_intersect << ',' << n_union << ',' << n_overlap << ',';
            fs_out << start << ',' << end << ',' << code_start << ',' << code_end << endl;
        }
    }
    cout << "processed " << encode_data.size() << " reads" <<endl;
    fs_out.close();
}

void Assembler::correct_reads(string encode_file, string align_file, string cmpreads_diff_file, string out_file)
{
    if (this->cand_size > 32)
        throw runtime_error("cand_size should not exceed 32.");
    
    // load encode and alignment files
    cout << "pileup encode_file" << endl;
    call_pileup_var(encode_file);
    cout << "pileup align_file" << endl;
    call_pileup_reads(align_file);
    
    // generate template
    counter = 1;
    temp_var = vector<uint64_t>(this->n_reads, 0);
    temp_read = vector<uint64_t>(this->n_reads, 0);
    
    // open input and output files
    FILE *p_infile = fopen(cmpreads_diff_file.c_str(), "rb");
    if (p_infile == NULL)
        runtime_error("fail to open cmpreads_diff_file");

    
    ofstream fs_outfile; open_outfile(fs_outfile, out_file);
    
    // scan cmpreads_diff_file
    CmpreadsDiffRead cur_cmpread(-1);
    int64_t n_lines = 0;
    while(true){
        // read line
        int cand_loci_size;
        int cand_loci_diff_size;
        int read_id;
        int start;
        int end;
        fread(&read_id, sizeof(int), 1, p_infile);
        fread(&start, sizeof(int), 1, p_infile);
        fread(&end, sizeof(int), 1, p_infile);
        
        fread(&cand_loci_size, sizeof(int), 1, p_infile);
        vector<int> cand_loci(cand_loci_size,-1);
        fread(&cand_loci[0], sizeof(int), cand_loci_size, p_infile);
        
        fread(&cand_loci_diff_size, sizeof(int), 1, p_infile);
        vector<int> cand_loci_diff(cand_loci_diff_size,-1);
        fread(&cand_loci_diff[0], sizeof(int), cand_loci_diff_size, p_infile);
        
        if (feof(p_infile)){
            // correct the currect read
            this->correct_reads_core(cur_cmpread);
            
            // output results
            this->print_correct_reads_raw(cur_cmpread, fs_outfile);
            
            break;
        }
        
        ++n_lines;
        if (n_lines % 1000 == 0)
            cout << "# of lines: " << n_lines << endl;
        
        // load CmpreadsDiff
        CmpreadsDiff cur_cmp(cand_loci, cand_loci_diff);
        cur_cmp.start = start;
        cur_cmp.end = end;
        
        // if get into a new read, correct the currect read
        if (read_id != cur_cmpread.read_id){
            if (cur_cmpread.read_id != -1){
                // correct the currect read
                this->correct_reads_core(cur_cmpread);
                
                // output results
                this->print_correct_reads_raw(cur_cmpread, fs_outfile);
            }
            // update read_id and clean cur_cmpread
            cur_cmpread.read_id = read_id;
            cur_cmpread.cmpreads_diff.clear();
            cur_cmpread.encode_corrected.clear();
        }
        
        // update cur_cmpread
        cur_cmpread.cmpreads_diff.push_back(cur_cmp);

    }
    cout << "# of lines: " << n_lines << endl;
    fclose(p_infile);
    //fclose(p_outfile);
    fs_outfile.close();
    
    // clean template and counter
    counter = 0;
    temp_var.clear();
    temp_read.clear();
    
}

void Assembler::correct_reads_core(CmpreadsDiffRead &cmpread)
{
    // scan each cmpreads_diff to test
    for (int i=0; i<cmpread.cmpreads_diff.size(); ++i){
        // test common variants
        for (int j = 0; j < (int)cmpread.cmpreads_diff[i].cand_loci.size(); ++j){
            // get focal locus
            int focal_locus = cmpread.cmpreads_diff[i].cand_loci[j];
            
            // get loci_set
            int win_start = j - (cand_size-1) > 0 ? j - (cand_size-1) : 0;
            for (int k = win_start; k <= j; ++k){
                // get loci_set
                int cur_end = k + (cand_size-1) < (int)cmpread.cmpreads_diff[i].cand_loci.size() - 1 ? k + (cand_size-1) : (int)cmpread.cmpreads_diff[i].cand_loci.size() - 1;
                if (k == cur_end) continue;
                vector<int> loci_set;
                for (int t = k; t <= cur_end; ++t){
                    if (t != j)
                        loci_set.push_back(cmpread.cmpreads_diff[i].cand_loci[t]);
                }
                
                // test conditional probability of focal locus given loci_set
                double logLR, condprob;
                int n_y_xp, n_xp;
                this->test_locus(focal_locus, loci_set, logLR, condprob, n_y_xp, n_xp);
                
                // record results
                if (condprob > cmpread.cmpreads_diff[i].condprob[j] &&
                    n_y_xp >= this->min_count){
                    cmpread.cmpreads_diff[i].condprob[j] = condprob;
                }
                
                if (cur_end == (int)cmpread.cmpreads_diff[i].cand_loci.size() - 1)
                    break;
            }
            
        }
        
        // test different variants
        for (int j = 0; j < (int)cmpread.cmpreads_diff[i].cand_loci_diff.size(); ++j){
            // get focal locus
            int focal_locus = cmpread.cmpreads_diff[i].cand_loci_diff[j];
            
            // get loci_set
            for (int k = 0; k < (int)cmpread.cmpreads_diff[i].cand_loci.size(); ++k){
                int cur_end = k + (cand_size-1) < (int)cmpread.cmpreads_diff[i].cand_loci.size() - 1 ? k + (cand_size-1) : (int)cmpread.cmpreads_diff[i].cand_loci.size() - 1;
                if (k == cur_end) continue;
                
                vector<int> loci_set;
                for (int t = k; t <= cur_end; ++t){
                    if (t != j)
                        loci_set.push_back(cmpread.cmpreads_diff[i].cand_loci[t]);
                }
                
                // test conditional probability of focal locus given loci_set
                double logLR, condprob;
                int n_y_xp, n_xp;
                this->test_locus(focal_locus, loci_set, logLR, condprob, n_y_xp, n_xp);

                // record results
                if (condprob > cmpread.cmpreads_diff[i].condprob_diff[j] &&
                    n_y_xp >= this->min_count){
                    cmpread.cmpreads_diff[i].condprob_diff[j] = condprob;
                }
                if (cur_end == (int)cmpread.cmpreads_diff[i].cand_loci.size() - 1)
                    break;
            }
        }
        
    }
    
}

void Assembler::test_locus(int focal_locus, const vector<int> &loci_set, double &logLR, double &condprob, int &n_y_xp, int &n_xp)
{
    // Let's fill in temp_vec_var and temp_vec_read
    // by reponse y. pu_var is 4 times larger than pu_read because of binary coding so
    // we have to devide 4 to access pu_read
    int y_locus = focal_locus;
    int y_read_locus = int (y_locus / 4);
    
    // fill in temp_vec_var and temp_vec_read by response y
    for (int j = 0; j < pu_var[y_locus].size(); j++)
        temp_var[pu_var[y_locus][j]] = counter;
    
    for (int j = 0; j < pu_read[y_read_locus].size(); j++)
        temp_read[pu_read[y_read_locus][j]] = counter;
    
    ++counter;
    
    if (counter >= numeric_limits<uint64_t>::max()-1)
        throw runtime_error("counter exceeds maximal int64_t");
    
    // calculate conditional probability
    for (int i = 0; i < (int)loci_set.size(); ++i){
        int cur_locus = loci_set[i];
        if (cur_locus == y_locus)
            continue;
        n_y_xp = 0; n_xp = 0;
        for (int k = 0; k < pu_var[cur_locus].size(); k++){
            if (temp_var[ pu_var[cur_locus][k] ] == counter - 1){
                temp_var[ pu_var[cur_locus][k] ] = counter;
                ++n_y_xp;
            }
            if (temp_read[ pu_var[cur_locus][k] ] == counter - 1){
                temp_read[ pu_var[cur_locus][k] ] = counter;
                ++n_xp;
            }
        }
        ++counter;
        if (counter >= numeric_limits<uint64_t>::max()-1)
            throw runtime_error("counter exceeds maximal int64_t");
    }
    
    if (n_xp == 0)
        condprob = 0;
    else
        condprob = (double)n_y_xp / n_xp;
    
}

void Assembler::run(string encode_file, string align_file, string out_file)
{
    

    
}


void Assembler::print_correct_reads_raw(const CmpreadsDiffRead &cmpread, ofstream &fs_outfile)
{
    for (int i = 0; i < (int)cmpread.cmpreads_diff.size(); ++i){
        fs_outfile << cmpread.read_id << '\t' << cmpread.cmpreads_diff[i].start << '\t' << cmpread.cmpreads_diff[i].end << '\t';
        fs_outfile << cmpread.cmpreads_diff[i].cand_loci << '\t' << cmpread.cmpreads_diff[i].cand_loci_diff << '\t';
        fs_outfile << cmpread.cmpreads_diff[i].condprob << '\t' << cmpread.cmpreads_diff[i].condprob_diff << endl;
    }
}


