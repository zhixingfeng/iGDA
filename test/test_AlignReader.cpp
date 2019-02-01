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

TEST_CASE("Test AlignReaderSam read()", "[hide]"){
    AlignReaderSam AlignReaderSam_obj;
    AlignReader *p_align = &AlignReaderSam_obj;
    p_align->getref("../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.ref.fa");
    p_align->open("../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.sam");
    Align align;
    ofstream fs_outfile;
    open_outfile(fs_outfile, "../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.fromsam");
    while(p_align->readline(align)){
        fs_outfile << align.qName << ' ' << align.tStart << ' ' << align.tEnd + 1 << ' ' << align.tStrand << ' ';
        fs_outfile << align.mapQV << ' ' << align.qAlignedSeq << ' ' << align.matchPattern <<  ' ' <<align.tAlignedSeq << endl;
    }
    fs_outfile.close();
    
    string cmd = "cut -f 1,9,10,11,17,18,19,20 -d ' ' ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5 > ../results/encode_from_sam/tmp.txt";
    system(cmd.c_str());
    
    cmd = "sort ../results/encode_from_sam/tmp.txt > ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.sorted";
    system(cmd.c_str());
    
    cmd = "sort ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.fromsam > ../results/encode_from_sam/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.toref.m5.fromsam.sorted";
    system(cmd.c_str());
    //p_align->readline(align);
    //cout << align.qAlignedSeq << endl;
    //cout << align.matchPattern << endl;
    //cout << align.tAlignedSeq << endl;
}

TEST_CASE("AlignReaderSam::samtom5()", "[hide]"){
    string samfile = "../results/NCTC3000/ERS950465_head.sam";
    string reffile = "../results/NCTC3000/ERS950465.fna";
    string m5file = "../results/NCTC3000/ERS950465_head.m5";
    AlignReaderSam alignreadersam;
    alignreadersam.samtom5(samfile, reffile, m5file);
    //alignreadersam.getref(reffile);
    
}



