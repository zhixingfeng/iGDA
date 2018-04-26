//
//  test_recode.cpp
//  iGDA
//
//  Created by Zhixing Feng on 2018/4/13.
//  Copyright © 2018年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/io.h"
#include "../src/misc/basic.h"
#include "../src/modules/assemble/assembler.h"


TEST_CASE("test assembler::correct_reads() (no recoding)", "[hide]")
{
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.5000";
    string cmpreads_diff_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.cmpreads_diff.5000";
    string out_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000.corrected";
    
    Assembler assembler;
    assembler.correct_reads(encode_file, align_file, cmpreads_diff_file, out_file, false);
}


TEST_CASE("test assembler::correct_reads() (recoding)", "[hide]")
{
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.5000";
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000";
    string cmpreads_diff_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.cmpreads_diff.5000";
    string out_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected";
    
    Assembler assembler(5, 20, 10,
                        0.2, 0.75);
    vector<int> read_ids = assembler.correct_reads(encode_file, align_file, cmpreads_diff_file, out_file, true);
    select_lines(read_ids, align_file, out_file + ".m5");

}


TEST_CASE("test assembler::olc() (recoding)", "[hide]")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected";
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected.m5.trimmed";
    string olc_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected.olc";
    
    Assembler assembler;
    assembler.olc(encode_file, align_file, olc_file, "");
}

TEST_CASE("test assembler::olc() (recoding, ncreads)")
{
    string encode_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected";
    string align_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected.m5.trimmed";
    string olc_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected.olc.nc";
    string follower_file = "../results/realign/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.recode.5000.corrected.non_contained.check_follower";
    
    Assembler assembler;
    assembler.olc(encode_file, align_file, olc_file, follower_file);
}












