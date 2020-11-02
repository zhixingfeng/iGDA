//
//  dforestsnvstxxl.cpp
//  iGDA
//
//  Created by Zhixing Feng on 10/31/20.
//  Copyright © 2020 Zhixing Feng. All rights reserved.
//

#include "dforestsnvstxxl.h"

// stxxl_vv_int
void stxxl_vv_int::pileup_encode(string encode_file)
{
    // clear data
    this->dat_vec.clear(); vector<stxxl_v_int>().swap(this->dat_vec);
    this->dat.clear(); stxxl::vector<int>().swap(this->dat);
    
    // get size pileup
    ifstream p_encode_file; open_infile(p_encode_file, encode_file);
    size_t pu_size = 0;
    while (true) {
        string buf; getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (auto i = 0; i < buf_vec.size(); ++i)
            if (buf_vec[i] + 1 > pu_size) pu_size = buf_vec[i] + 1;
    }
    p_encode_file.close();
    if (pu_size == 0) return;
    
    // get size of each locus
    this->dat_vec = vector<stxxl_v_int>(pu_size, stxxl_v_int(&this->dat));
    size_t total_size = 0;
    open_infile(p_encode_file, encode_file);
    while (true) {
        string buf; getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (auto i = 0; i < buf_vec.size(); ++i)
            ++this->dat_vec[buf_vec[i]].size();
        total_size += buf_vec.size();
    }
    p_encode_file.close();
    if (total_size == 0) return;
    this->dat.resize(total_size);
    
    // construct shift of stxxl_v_int
    this->dat_vec[0].shift = 0;
    for (auto i = 0; i < this->dat_vec.size() - 1; ++i){
        this->dat_vec[i+1].shift = this->dat_vec[i].shift + (int) this->dat_vec[i].size();
    }
    
    // pileup on hard drive
    vector<int> pu_shift(pu_size, 0);
    open_infile(p_encode_file, encode_file);
    int read_id = 0;
    while (true) {
        string buf; getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (auto i = 0; i < buf_vec.size(); ++i){
            this->dat[ this->dat_vec[buf_vec[i]].shift + pu_shift[buf_vec[i]] ] = read_id;
            ++pu_shift[buf_vec[i]];
        }
        ++read_id;
        if ( read_id >= std::numeric_limits<int>::max() - 1 )
            throw runtime_error("stxxl_vv_int::pileup_encode() exceed maximal number of reads " + to_string(std::numeric_limits<int>::max() - 1) + "." );
    }
    p_encode_file.close();
}

void stxxl_vv_int::print_dat_vec(string outfile)
{
    ofstream fs_outfile;
    open_outfile(fs_outfile, outfile);
    for (auto i = 0; i < this->dat_vec.size(); ++i){
        if (this->dat_vec[i].size() > 0){
            for (auto j = 0; j < this->dat_vec[i].size(); ++j)
                fs_outfile << this->dat_vec[i][j] << "\t";
            fs_outfile << endl;
        }
    }
    fs_outfile.close();
}
