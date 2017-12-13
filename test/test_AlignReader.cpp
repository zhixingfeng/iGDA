//
//  test_AlignReader.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/4/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/modules.h"

TEST_CASE("Test AlignReaderM5", "[hide]"){
    AlignReaderM5 AlignReaderM5_obj;
    AlignReader *p_align = &AlignReaderM5_obj;
    Align align;
    
    p_align->open("../data/MSSA_61_forward.m5");
    
    while (p_align->readline(align)){
        
    }
    p_align->close();
}

TEST_CASE("Test AlignReaderM5 read()", "[hide]"){
    AlignReaderM5 AlignReaderM5_obj;
    AlignReader *p_align = &AlignReaderM5_obj;
    
    // read m5 file
    stxxl::vector<Align> align_vec;
    p_align->read("../data/MSSA_61_forward.m5", align_vec);
    
    // write m5 file
    ofstream fs_outfile; open_outfile(fs_outfile, "../data/MSSA_61_forward.m5.qAlignedSeq.reconstruct");
    for (int i=0; i<(int)align_vec.size(); ++i){
        fs_outfile << align_vec[i].qAlignedSeq << endl;
    }
    fs_outfile.close();
}
