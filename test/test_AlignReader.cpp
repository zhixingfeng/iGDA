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

#include <seqan/sequence.h>
#include <seqan/bam_io.h>

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
    vector<Align> align_vec;
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

TEST_CASE("AlignReaderSam::samtom5qv()", "[hide]"){
    string sam_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/pileup/SAMEA4916110_split/NZ_UWXV01000001.1_forward.sam";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/pileup/GCF_900608245.1_kpneu039_genomic_NZ_UWXV01000001.1.fna";
    string m5qv_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/ont_kp/pileup/SAMEA4916110_split/NZ_UWXV01000001.1_forward.m5qv";
    AlignReaderSam alignreadersam;
    alignreadersam.samtom5qv(sam_file, ref_file, m5qv_file);
    
}

TEST_CASE("AlignReaderSam::bamchrrange()", "[hide]")
{
    //string sam_file = "/Volumes/zxfeng_no_1/work/igda_memory/pacbio_meta/data/merged_sorted_chr_qv8.bam";
    //string ref_file = "/Volumes/zxfeng_no_1/work/igda_memory/pacbio_meta/reference/Borrelia_burgdorferi_B31_chr.fna";
    string sam_file = "/Users/zhixingfeng/Dropbox/work/iGDA/paper/submission/nature_communication_revision/analysis/memory_reduce/pacbio_ecoli/data/merged_realign_trim.bam";
    string ref_file = "/Users/zhixingfeng/Dropbox/work/iGDA/paper/submission/nature_communication_revision/analysis/memory_reduce/pacbio_ecoli/reference/ecoli_K12_MG1655.fasta";
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, "/Users/zhixingfeng/Dropbox/work/iGDA/paper/submission/nature_communication_revision/analysis/memory_reduce/pacbio_ecoli/results/igda/detect/qv0/realign.m5");
    AlignReaderSam alignreadersam;
    alignreadersam.getchrrange(sam_file, ref_file, sam_file + ".chrrange");
    
}
