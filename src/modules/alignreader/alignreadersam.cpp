//
//  alignreadersam.cpp
//  iGDA
//
//  Created by Zhixing Feng on 6/29/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "alignreadersam.h"


bool AlignReaderSam::open(string filename) {
    p_file.open(filename);
    if (!p_file.is_open())
        throw runtime_error("fail to open " + filename);
    return true;
}

bool AlignReaderSam::readline(Align &align) {
    
    return true;
}

bool AlignReaderSam::close() {
    p_file.close();
    return true;
}


bool AlignReaderSam::read(string filename, stxxl::vector<Align> &align_vec)
{
        
    return true;
}
