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


TEST_CASE("test stxxl", "[hide]")
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


TEST_CASE("test vector of vector")
{
    stxxl::vector<vector<int> > cmpreads;
    // read
    ifstream fs_infile; open_infile(fs_infile, "../results/B_10_cons_cmpreads.txt");
    while(true){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof())
            break;
        cmpreads.push_back(split_int(buf, ','));
    }
    fs_infile.close();
    
    // write
    clock_t time_begin = clock();
    ofstream fs_outfile; open_outfile(fs_outfile, "../results/B_10_cons_cmpreads.txt.stxxl");
    for (int i=0; i<(int)cmpreads.size(); ++i)
        fs_outfile << cmpreads[i] << endl;
    fs_outfile.close();
    clock_t time_end = clock();
    cout << "time elapse: " << double(time_end - time_begin) / CLOCKS_PER_SEC << endl;
}




