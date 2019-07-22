//
//  test_rsm.cpp
//  iGDA
//
//  Created by Zhixing Feng on 7/22/19.
//  Copyright © 2019 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../src/modules/rsm/rsmsnv.h"
#include <ctime>


TEST_CASE("test rsmsnv", "[hide]")
{
    string encode_file = "../data/nanopore_kp/igda_pipe_maskqv/merged_maskqv.encode";
    string encode_ref_file = "../data/nanopore_kp/igda_pipe_maskqv/merged_maskqv.encode.ref";
    string m5_file = "../data/nanopore_kp/igda_pipe_maskqv/merged_maskqv.m5";
    string cmpreads_file = "../data/nanopore_kp/igda_pipe_maskqv/merged_maskqv.cmpreads";
    string ref_file = "../data/nanopore_kp/igda_pipe_maskqv/k_pneumoniae_hs11286_chr.fasta";
    string out_file = "../data/nanopore_kp/igda_pipe_maskqv/merged_maskqv.rsm";
    
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    RSMsnv rsmsnv(&alignreader, &aligncoder);
    rsmsnv.load_homo_blocks(ref_file);
    rsmsnv.run(encode_file, encode_ref_file, m5_file, cmpreads_file, out_file, 25, 10000);
    
}
