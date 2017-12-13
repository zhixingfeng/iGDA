//
//  alignreader.h
//  iGDA
//
//  Created by Zhixing Feng on 8/4/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef __iGDA__alignreader__
#define __iGDA__alignreader__

#include <headers.h>
#include <stxxl.h>
#include "../../misc/seqopt.h"
// position in m5 format is 1-based and (], and convert to 0-based [] when reading
struct Align
{
    // original features in m5 format
    string qName;
    int qLength;
    int qStart;
    int qEnd;
    char qStrand;
    string tName;
    int tLength;
    int tStart;
    int tEnd;
    char tStrand;
    int score;
    int numMatch;
    int numMismatch;
    int numIns;
    int numDel;
    int mapQV;
    string qAlignedSeq;
    string matchPattern;
    string tAlignedSeq;
    
    // additional feature (query seq and ref seq without not "-". Used for realignment)
    string qSeq;
    string tSeq;
};

class AlignReader
{
public:
    AlignReader(){}
    virtual ~AlignReader(){}
    
    /*-----------virtual functions------------*/
    // open align file
    virtual bool open(string filename)=0;
    
    // read a line of align file
    virtual bool readline(Align &align)=0;
    
    // close align file
    virtual bool close()=0;
    
    // read all alignment and store it into stxxl vector
    virtual bool read(string filename, stxxl::vector<Align> &align_vec)=0;
        
    
protected:
    string filename;
};


#endif /* defined(__iGDA__alignreader__) */
