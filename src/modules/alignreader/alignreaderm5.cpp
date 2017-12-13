//
//  alignreaderm5.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/4/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "alignreaderm5.h"

bool AlignReaderM5::open(string filename) {
    p_file.open(filename);
    if (!p_file.is_open())
        throw runtime_error("fail to open " + filename);
    return true;
}

bool AlignReaderM5::readline(Align &align) {
    string buf;
    getline(p_file, buf);
    if (p_file.eof())
        return false;
    vector<string> buf_vec = split(buf, ' ');
    if ((int) buf_vec.size() != 19)
        throw runtime_error("incorrect format in " + filename);
    
    align.qName = buf_vec[0];
    align.qLength = stod(buf_vec[1]);
    align.qStart = stod(buf_vec[2]);    // range in m5 is (], use 0-based position
    align.qEnd = stod(buf_vec[3]) - 1;
    align.qStrand = buf_vec[4][0];
    align.tName = buf_vec[5];
    align.tLength = stod(buf_vec[6]);
    align.tStart = stod(buf_vec[7]);    // range in m5 is (], use 0-based position
    align.tEnd = stod(buf_vec[8]) - 1;
    align.tStrand = buf_vec[9][0];
    align.score = stod(buf_vec[10]);
    align.numMatch = stod(buf_vec[11]);
    align.numMismatch = stod(buf_vec[12]);
    align.numIns = stod(buf_vec[13]);
    align.numDel = stod(buf_vec[14]);
    align.mapQV = stod(buf_vec[15]);
    align.qAlignedSeq = buf_vec[16];
    align.matchPattern = buf_vec[17];
    align.tAlignedSeq = buf_vec[18];
    
    return true;
}

bool AlignReaderM5::close() {
    p_file.close();
    return true;
}


bool AlignReaderM5::read(string filename, stxxl::vector<Align> &align_vec)
{
    Align align;
    
    this->open(filename);
    while (this->readline(align)){
        align.qSeq = rm_indel_from_seq(align.qAlignedSeq);
        align.tSeq = rm_indel_from_seq(align.tAlignedSeq);
        align_vec.push_back(align);
    }
    this->close();

    return true;
}
