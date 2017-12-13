//
//  test_stxxl.cpp
//  iGDA
//
//  Created by Zhixing Feng on 12/13/17.
//  Copyright (c) 2017 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/io.h"
#include "../src/misc/basic.h"
#include "../src/modules/assemble/assembler.h"
#include "../src/modules/alignment/alignment.h"
#include "../tools/tools.h"
#include <ctime>


TEST_CASE("test stxxl")
{
    typedef stxxl::VECTOR_GENERATOR<int>::result vector;
    vector my_vector;
    for (int i = 0; i < 1024 * 1024; i++){
        my_vector.push_back(i + 2);
    }
    std::cout << my_vector[99] << std::endl;
    my_vector[100] = 0;
    while (!my_vector.empty())
    {
        my_vector.pop_back();
    }

}