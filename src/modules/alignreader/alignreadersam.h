//
//  alignreadersam.h
//  iGDA
//
//  Created by Zhixing Feng on 6/29/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__alignreadersam__
#define __iGDA__alignreadersam__

#include "alignreader.h"

class AlignReaderSam : public AlignReader
{
public:
    AlignReaderSam(){filename="";}
    virtual ~AlignReaderSam(){}
    
    // open align file
    bool open(string filename);
    
    // read a line of align file
    bool readline(Align &align);
    
    // close align file
    bool close();
    
    // read all alignment and store it into stxxl vector
    bool read(string filename, stxxl::vector<Align> &align_vec);
    
protected:
    string filename;
    ifstream p_file;
};

#endif /* defined(__iGDA__alignreadersam__) */
