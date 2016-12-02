//
//  io.h
//  iGDA
//
//  Created by Zhixing Feng on 16/10/9.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef io_h
#define io_h
#include "../../include/headers.h"
#include "../modules/alignreader/alignreaderm5.h"

typedef pair<int,int> ReadRange;

inline bool loadencodedata(vector<vector<int> > &encode_data, string encode_file)
{
    ifstream p_encode_file; open_infile(p_encode_file, encode_file);
    while (true) {
        string buf;
        getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        encode_data.push_back(split_int(buf, '\t'));
    }
    p_encode_file.close();
    return true;
}

// format: m=m5
inline bool loadreadsrange(vector<ReadRange> &reads_range, string align_file, char format='m')
{
    // setup alignreader
    AlignReader *p_alignreader;
    AlignReaderM5 alignreaderm5;
    switch(format){
        case 'm':
            p_alignreader = &alignreaderm5;
            break;
        default:
            throw runtime_error("loadreadsranges: unsupported format.");
    }
    
    // load alignment data
    Align align;
    p_alignreader->open(align_file);
    while(p_alignreader->readline(align))
        reads_range.push_back(ReadRange(align.tStart, align.tEnd));
    
    p_alignreader->close();
    
    return true;
}


#endif /* io_h */
