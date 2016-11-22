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