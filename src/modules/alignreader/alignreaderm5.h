//
//  alignreaderm5.h
//  iGDA
//
//  Created by Zhixing Feng on 8/4/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__alignreaderm5__
#define __iGDA__alignreaderm5__

#include "alignreader.h"

class AlignReaderM5 : public AlignReader
{
public:
    AlignReaderM5(){filename="";}
    virtual ~AlignReaderM5(){}
    
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


#endif /* defined(__iGDA__alignreaderm5__) */
