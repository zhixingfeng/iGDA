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
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>


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
    
    // get reference
    bool getref(string filename);
    
    // convert sam to m5
    bool samtom5(string sam_file, string ref_file, string m5_file);
    
protected:
    string filename;
    seqan::BamFileIn bamFileIn;
    seqan::StringSet<seqan::CharString> ref_ids;
    seqan::StringSet<seqan::Dna5String> ref_seqs;

};

#endif /* defined(__iGDA__alignreadersam__) */

