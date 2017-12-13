//
//  test_AlignCoder.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/misc.h"
#include "../src/modules/modules.h"
#include "../src/modules/alignment/alignment.h"
TEST_CASE("test AlignCoderSNV::encode", "[hide]")
{
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->encode("../data/MSSA_61_forward.m5", "../results/MSSA_61_forward_encode_snv.txt");
}

TEST_CASE("test AlignCoderSNV::encode negative strand","[hide]")
{
    
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;
    AlignReaderM5 alignreaderm5;
    p_aligncoder->setAlignReader(&alignreaderm5);
    p_aligncoder->encode("../data/B_10_cons.m5", "../results/B_10_cons.encode");
}


TEST_CASE("test AlignCoderSNV::decode", "[hide]")
{
    AlignCoderSNV aligncodersnv;
    AlignCoder *p_aligncoder = &aligncodersnv;

    REQUIRE(p_aligncoder->decode(2).first==1); REQUIRE(p_aligncoder->decode(2).second=='C');
    REQUIRE(p_aligncoder->decode(3).first==1); REQUIRE(p_aligncoder->decode(3).second=='G');
    REQUIRE(p_aligncoder->decode(4).first==1); REQUIRE(p_aligncoder->decode(4).second=='T');
    REQUIRE(p_aligncoder->decode(5).first==2); REQUIRE(p_aligncoder->decode(5).second=='A');
    REQUIRE(p_aligncoder->decode(6).first==2); REQUIRE(p_aligncoder->decode(6).second=='C');
    REQUIRE(p_aligncoder->decode(7).first==2); REQUIRE(p_aligncoder->decode(7).second=='G');
    REQUIRE(p_aligncoder->decode(8).first==2); REQUIRE(p_aligncoder->decode(8).second=='T');
    
    REQUIRE(p_aligncoder->decode(3894).first==974); REQUIRE(p_aligncoder->decode(3894).second=='C');
    REQUIRE(p_aligncoder->decode(372).first==93);REQUIRE(p_aligncoder->decode(372).second=='T');
    
    
    
}



TEST_CASE("test AlignCoderSNV::encode (ssw)", "[hide]")
{
    string read = "AGGGCGGGGTTGTCCCAAGCCCGGCGCGGGGCGAAGCGGGTCGCCGTGGCTGCTGCAGCGCCTTCAGCCGGGCTCGGTACATCGCGTGTGTGCTCTCTCGTGTGTGAGTTGCGACGCAGCGGATGTCTGGATCGGACCGCCGCGGAATCGGCGCGAAAGCAGTACGAAGATCACGGCCAGCATCCAGCAGGTCCGCGTCGGGCCAGCGCACGCCCTCGTGCAGGCAGGCCCGCTTCGTTCAACAGCGGCAGGATCAGGCGTCCTCCAGATCGGCAGGGGCCTGGTAGTCACTGGGGCCACGTCCGGCTTCTTCGGGCTTCCGTTCTCCCAGGGTGTAGTGCCCTGCCATCCTTCTTGCCGCGCTTGCGCGGCTTTCAACCGTCTGCAGGGCGGCCGGGATCTGCAGGCCCAGGAACGGCGCAGTTCGCGGCCGCGCCGGGCGGCACGTCCAGGCCACCGTGGTCGATCAGCTCGATCGGGCCCATCGGCATGGCCGAACTTGCACCGCAGCCCTGGGTCGATCACCGGGCCCGGAATGGCCTCGGCATAGGCGGTGGCGCTTCCAGCATGTACCGGGAACAGCACGCGGTGCCAGGAAGGCCGGCTGCCGGGCCAACCGCACGGGCAACGTTGCCCAGGCCTTGCAAACGCGGCCAGGCGGGCGCTCGGTTTCCGGCGCCATGCCGTCTGCGTGGATGATCTCCCCCCCCCACCAGCGGCACTGCGCGAACCGGGTGAATTAGTGCAAGGCCAGCGAACTGCGCCGGGCGCTGGATGGTGGTCGCGCAGTTCCACCAACGGAATCGAGGAGGATGTGGTGGTCGCAGCGCGTCAGCTTCATCTTCGGTTTTCCAGCGTCTGGTACAGCGCACGCTTGCCTCAGGGTTCTCGATGAATCGCTTCGATGCACCAGGTCGGCCGTCGGCCACACCGTGGCCTTCCAGGTCGGCGCGGCAGTGCGCGCGGCCACGGCACGGGCGCTTGCTTTCATCGCGCACCTTCTTGGCGAACAGCGCCTGCGCGCGCTCCCCATGCCGGGTCAGAGCGCTGCGTCGCGGTCCTGCAGGTCACCTGCGAAGCCCTTGTAGGCGGCCCCATGCGGCGCGATGTCACCGCCACATCACCGCCGGCACCGGACCACGGTGCACCGTGGCGGATGGCCGGGAATCACCGCCGCCCAGGCCCTTCAGGCGTTCTGGTCAGGAAGAAGATGCGGATCAGGTTGCGCGCGGTCGGCGTGCTCGCCAGTTCACCCACCGCACGGCGTTCGGCGTCCAGGCGGCCTGGATCGGCTTGCCGCCGCTGCGCTGCCAGGGTGCTGATCAGCGCATACGGCGCCGGATACTGGTCCTTCTTCGCCTTGCGCGCGAACCTGCTTGACCATCTGCGGGGCCAGGCAGCGTGCGCGCCGAGCCAGGTGTTGGTGGCCCACGCGGTGGCGCGCTGCGTTGAACGGGGGTGGGTGCCGGCAGGCAGCGCCACGGCGTGTCGAGCACCCACTGCCGGTGCCACCACCTTGTCGACCCAGGCCGATGCCACGCGCGGCCGAGGCCGGACAGGGTACGGGCCGGTCAACACATTCATCGCGTGCGGGAGCGCCACCAGGCTTGGGCAGGCGGGCGCTGCCGCCCCAGCCGGGGAAGATGCCCAGCTGGGTTCCGGCAGGCCGGATGGCGGGTGCTGCTGTCATTGGACGCCACGGCGTACGGCAGGCCAGCTCGCAGGCTCGGTGCCGCCGCCCAGGCAGTGGCCGTGGAATGGCCGCCACCGCTAGGGCAGGCAGCTCGGCCAGCGTTCTGGTACGTGCCTGGCGCGACGGATCGCGTCGGTTGACGGTGCCGCGCGATCGAATTCCTGGAACCCTTCAGTCGGCCCGGCAATGAAAGCGCCTTCTTCAACGAACTGGATGACCACGTCCCTTTGGGCGGGTCCAGAGGCAATGCGCTGCAAAGCAGGTCGCCCAGTTCCAGCAGCACGTCCTTGCCGACATCGCATTTGACGCTGCTGTTCCCGGCGATCCAGGGAGAGAACCAGCCACGCGCGTCGTCTGCGGATCTCGACGGTGCAGTGGCTGAAGCGGAGACCATCGAGCCTGAGAGCATGGACGTTCCTCCGGTTGGTGTATGGAGAAGCGGCTATTGATCCAGAGGTTGTTTCGCCCGCTCAGAATCCCAATAGCAGGGGGGTGTGTCGTTCCCCGGACCATGGCCGTGGTACCGGAACGCGCTAACGGAGACCGGTGTTAAAAGTGGTGCGGCGCCAGTGATCGGCGACTGACTTTTTCCCGGGGATTGGCAGTCTCATTGGCCTACGTCGGTCTGGGCTGACAGAGTGCGGGCCCCTCATGGCCGATGTTGGATACACCGTCAGGAAGTGGATTGGGGAGCTGGGGTCCGGCGCGTGCAGCAACGGCGAGGGAGTGCCGCGGTTCGATGCCCTGGTGGCGCAAGTCCAGCACCGGTCGGTCGTGGCCCTGGTCGGTCGTCTACATCCGCCGACTGAGCGAAGTGTCAGGACGTCGCCAGGACACTTCATCCGTGCTACCGCGCGATCGGGAGTTTCCGTGCGATGCCCAGTCTCAAACCTGGTTGCATGCGAATCCGTGAGATCCGCCAAAGACTAGCCTGGGTTCACAGCAGGATCGACGGGCCGCCTACCGGATGACGCGACTCGAGGCGATCGGCGATGGCCGAAACAGTTGGACAGGCGGTACCCCGCCTGCGCGACACCGATACGCCCCGGAACGCGAGTTGATGCGCCAGGAACTGGAACAGACGGTCATGAAGCCGTCAACGCGCATGCCGGAGAACTCCGGTCGGCGATCACCTCTGCGCGCGAGGTGGGAAGCCTGAGCTACGAGTATCTCGCGCAGAGAGATGGGGGTGCGGGCCGATCGGCACCGTGCGTTCACGGATCTTCCGGTCGCGCGGGCGATCGACACCGAACTCCGGCCGCTGTTGGAATCGCAAGCTGCCACCCGTGGAGAAGAGCGCGCGTATGACCAGTAACCCGGTTCACGAAGTCGCAGAACATCAATCGCCCGCCGGAGCAGCGCCTGGACCAGCGCCATCGCGAGCAGCTGTCGGCGCTGGTCGATGGCGAGCTCGGTTGCCGATGAAGGCGCTTCCTGCTACGCCGCATGAGCACGACCCCGAACTGGCCGGCTGCCAGGGAGCGCTGGCAGCTGCTGGGCGACGTGATGCGCGGGCAGCCCGGCGGCTGGCGCCGGCCGGCTCAGCGCCGCCGTGGCGCTGCCATCGCCGCCGAGCGGTGCCACCGTCCGAGCCGCGTCGCCAGGTGCGCCGCGGTGGCTGGCGGGCCTGGGGGGCGCGCCGCCCTGGCCGCGTCGGTGGCCGCCTCGCGCGTTCATCGGCAGCGTGAGAAACTGAAGGAGGCGGACGCCTGGTGGAGCCGCTGGCGCCCGCAGGTGATCGCCAGCCAGGCGGAACTGGCGCCGGCGCGGACCGCACCGCGCGGGCGGGGACCGAGCTCGGTGGATACGCGGATGGCCGGTGGTGCTGCGCCGGCCGTGGGCGGATTGGCCGCCTCCAGCCGCCGCCAGGACGCCCGCCGTGCCACGCCACGTACCAGAGTCCGGGCGTGCCCGCCGCAGCCGATGGGAACGGCCCCGCAACGGGCCATCGCCGCACCGGCGGCCGTTGACCGCCGGACCGTCCGGCCAACGGCCAACCGCAACCTGCCTTCGGCGAGTGGGGATTGCAGGCCGCGGCCGTGGCCGCGCTCCAGCCTGGCGCCCGCTGCCGGTGGCGCACTCAATGCCAGCTTCCCGGCCATGGCGGCGGCGCGGCGTTCTTACCGTTCGAGCGCGCCGTGCACGGGACGAGCTGCCGGTGCCGCGCGCGCCGCCGACGAGCCGCCTTTCACGCCCCTGCCGTCGTGCGGGCCGTTATTCCTCTTTCCCATCCGACCCGGAGGTGCCGTCGATGATCCCGACTTCGCACGCTAGCCATTGGCCTGCTGGCGCGCCTGAGCCGTTCGCCGTGTCGGCCTGCGCGCAGGCGGCCGGCGCCCGTCGCGCAGACCCAAGCCCGCGCCAGCCGCGCCGCGCGGCACCGGCGCAGCCTGGTCAGTGGCCTGCCCGACTTCCACCAACCTGGTCGAACAGGTCGGCCGGCTGTGTTCAAGCGTGCGATACCACCATCGTCCGCAGACAATCGCCAGCGCGCGCGGCGGCCGATGGGGCGACGATGACATGCCGGAGTTCTTCGCCGCTTCTTTGGTCCCACTTCCCGATGCCGGGGCAGGCCGGGTGGGCCGGACGCTGGCCCCGCATCAAGGGCCGTTGCATGGGCTCGGCTGTCATCATCTCGCCCGATGGCTATGTCGCTGACCAACTACCATGTGGTGGCGCGACGCCGTGCGACGTGGAGAGGTCAAGCTGGGCGACGATGGCCGTGAGTTCACCGCAAGGTGCGTGGGGCAGCGACCAGAGTTGACGACGTCGCGCTGCTGAAGATCGATGGAAGAACCTGCGACCGTGCGCGTGGGGGGTTCCAACACCCTTGAACCGGGCCAATGGGTCGTCGCGATCGGCTCGCCGTTCGGGCTCGACCATTCCGGTTACCGCGGCGTGGTCAGCGCCGTCGGCCGCAGCACGGGGCCCGACAGCGCTAGCGTGCCGTTCATCCAGACCGACGTGGCGGCAACCAGGCATTCCGGTGGCCCCTGCTGAATACCCGCGGCGAGGTGGTCGGCATCAACTCGCAGGATCTTCTCCCCTCCGGCGGGCTTCGATGGGCATCAGCTTCGCGATCCCGATCGATCTGGCGATGAGTGCGGTCGAGCAGATCAAGAGAGCGGCAAGTCACCGGGCCAGCTCGGCGCGGTGTTCGAGTCCGACGATTCGCTGAAGGCGCAGGGCCTGGGCCTGCCGGACAGCCGTGGCGCCCTGGATCACCAGATCGTGGCCCAGTGCGGCCCCAAGGCGGCGTGGCAGTGGGCGATTCATCCCGGCTCGGTCAACGGCGCCCGGTCAACAGGCTGGTCCGACTCTGCCGGCCGCATTGATTCGGTGCGAATGGCGCCGGGCGCAAGTGCAGACCTGGTGGTTGTACCGCGACGGGCAAAGCCGCGATGGACCTGCACGGCGCGACCCTGACCGACTGAGTGAAGACCGGCCAGCCCTGCCGCGGGCCGGCGGCCGGCGCCGAGGCGGCGCCGCAGACCGCGTCAATGCCCTGCTGGGGCTGGACTGTGACGCGACCTGAGCGCGGCCCAGCGCAAGCAAGTTCGGGCTCAACAGCAACGAAGGCGTGCGGACACCGCGTCAAGGCCAGCCGCTCGTGATGCCGGGCTGTCGCCGCGATGGTTGAGGTCCTGCAGTTGCCGTCGCGGCGGTCGGGCAGCGTGATGCCCTGAACCGGCGCTGTCCGCGTTCAAGGAAGGACGACGTGGTGATGCTGCTGGTCCGCACGGCCAACGGCAACAGCCGCCTCGGTGGGGTCCAAAGGCCGCCATAAGCCCGCCGGGGCCGCCCGCGCGCGTCTCCCCCCGAAGCCCCGTCCTGCGGGCTGTCCCGTTCCAGGGCGGGCACGATCGCGGGCCGTTTTCGCGCGGGGCGCCGGCATTGCGATAATCGACCGTTGACCCTGACGACGCCGCGCGCCGCCGCCACTATGTCTTCTGATTCGATGCGAACGGATCCGCAACTTCCCATCATCGGCCATGTCGACGCACGGCAAGTCCACGCTGGCGACCGCATCATCCAGCTTCGTGTGGTGGGCTGCAGGCCCGCGAGATGAAGGCAGGTGGCTCGACTCCAATCCGATCGAGCGCGAAGCTGGCATCACTTATCGAAGGCGAGCAGTCGGTGTCCGCTGCCCTACCTGGCCAAGGACGGGCAGACTACCACCTGAAGCTTCATCGACAGCCCCCGGGCGATTGTGCGACTTCTCGCTATGAAGTCAGCCGCTCGCTGGCCGCCTGCGAAGGCGCCGCTCTGGGTGGTCATGCGGCGCAGGGGCGTGGAAGGCCCAGTCGGTGGCCACTGCTACACCGCCGTGGAGCAGGGCCTGAAGTGCTGGGGTGATCAA";
    string ref = "AGGGCGGGGTTGTCCCAGCCCGGGCGCGGGGCGAAACGGTCGCCGTGGCGCTGCTGCAGCGCCTTCAGCCGCTCGACGATCGCATCGGCACCGACCGCGCGGATGTGCTGGATCGGACCGCCGCGGAACGGGCGAAACCGGTACCGAAGATCACGCCGGCATCCAGCAGGTCCGCATCGGCCACCACGCCCTCGTGCAGGCACGCCACCGCCTCGTTCAACAGCGGCAGGATCAGGCGGTCCTCCAGATCGGCCGGGGCCTGGTAATCGCTGGCCACGTCCGGCTTCTTCGGCTTGCCGTTTCCCAGGTGTAGATGCCCTGGCCGTCCTTCTTGCCGCGCTTGTCCGGTTCCACCGTCTGCAGTGCGGCCGGAATCTGCAGGCCCAGGAACGGCGCCAGTTCGCGGCCGACGCCGGCGGCCACGTCCAGGCCGACGGTGTCGATCAGTTCGATCGGGCCCATCGGCATGCCGAACTTCACCGCGGCCCTGTCGATCACCGGGCCCGGAATGCCCTCGGCATACGCGGTGGCGCTTCCAGCATGTACGGGAACAGCACGCGGTTGACCAGGAAGCCCGGGCTGCCGGCCACCGGCACCGGGAATTGCCCAGCGCCTTGCAGAACGCGGCCAGGCGGCGCTCGGTCTCCGGCGCCATGCCGTCGTGGTGGATGATCTCCACCAGCGGCATCTGCGCCACCGGGTTGAAGTAGTGCAGGCCAGCGAACTGCGCCGGGCGCTGGATGTGGTCGCGCAGTTCCACCAGCGGAATCGACGAGGTGTTGGTGGTCAGCAGCGCGTCCAGCTTCATCTTCGGTTCCAGCGTCTGGTACAGCGCACGCTTGGCTTCCGGGTTCTCGATGATCGCTTCGATCACCAGGTCGGCCTCGGCCACGCCGTTGCCTTCCAGGTCGGCACGCAGGCGCGCGGCCACGGCGGGCGCTTGCTCTCGTCGCGCACCTTCTTGGCGAACAGCGCCTGCGCGCGCTCCATGGCCGGGTCGATGAAGCGCTGCTCGCGGTCCTGCAGGGTCACCTCGAAGCCCTTGTAGGCGGCCCACGCGGCGATGTCACCGCCCATCACGCCGGCGCCGACCACGTGCACGTGGCGGATGCCGGAATCACCGCCGCCCAGGCCCTTCAGGCGTTCGGTCAGGAAGAAGATGCGGATCAGGTTGCGCGCGGTCGGCGTGCTGGCCAGCTTCACCACGCGCGGCGTTCGGCATCCAGGCGCGCCTGGATCGGCTTGCCGCCGCTGCGCTGCCAGGTGCTGATCAGCGCGTAGGCGCCGGATACTGGTCCTTCTTGCCTTGCGCGCGACCTGCTTGACCATCTGCGGTGCCAGCAGCGTGCGTGCCAGCCAGGTGTTGGTGGCCCAGGCCGTGGCGCGCTGCTTGAACGGGCGGGTGGTGCCGGACAGCGCCAGCGCCACGGCGGTGTCGAGCACCACGCCGGTGCCACCACCTTGTCGACCAGGCCGATGCCACGCGCGGCCGAGGCCGACAGGGTACGGCCGGTCAGCATCATGTCCATCGCTGCGGGGCGCCCACCAGCTGCGGCAGGCGCGCGCTGCCGCCCCAGCCGGGGAAGATGCCCAGCTGGGTTTCCGGCAGGCCGATGCGGGTGCTGCTGTCATTGGAGGCCACGCGGTAGCGGCAGGCCAGCGCCAGCTCGGTGCCGCCGCCCAGGCAGTGGCCGTGGATGGCCGCCACGGTGGGCAGGGCAGCTCGGCCAGCTTCTGGTAGTGGCCTGGCCGCGGCGGATGCGTCGTTGACGGTGCCGCGGCGGTCGAATTCCTGGAACTCCTTCAGGTCCGCACCGGCAATGAAGCCGGCCTTCTTCAGCGACTGGATCACCACGCCCTTGGGCGGGTCCAGTGCAATGCGCTCAAGCAGGTCGCCCAGTTCCAGCAGCACGTCCTGCGACATGCATTGACGCTGCTGTCCTGGCGATCCAGGGAGAGAACCACCACGCCGTCGTCGCGGATCTCGGGGTGCCAGTGGCTGAAGCGGAGACCATCGAAGCCTGAGAGCATGGACGTTCCATCCGGTTGTGTATGGAGAAGCGGCTATGATCCAGAGGTTGTTTCCCCCGCGTCAAATCCCAATAGCAGGGGAGTGGCGACCCCGGATCCATGGCCGTGGACCGGAACGCTAACGGAACGGTGGTTAAAAGCTGGTGCGGCGCCAGGTGATCGGCGACTGAACTTTTCCCGGGATTGGCAGTCTCATTGCCCGACGTCGGTTGGGCTGACAGGAGTGCGGCCCCTCATGGCCGATGTTGATACACCTCAGGAGCTGGATCTGGAACTGGTCCGGCGCGTGCAGCACGGCGAGAGCGCCGCGTTCGATGTCCTGGTGCGCAAGTACCAGCATCGGGTGTGGCCCTTGTCGGTCGCTACATCGCCGACTGGAGCGAATGTCAGGACGTCGCCCAGGACACTTTCATCCGTGCCTACCGCGCGATCGGAAGTTTCCGTGGCGATGCCCAGTTCTCAACCTGGTTGCATCGAATCGCCGTGAATACCGCCAAAAACTACCTGGCTTCACACAATCGACGGCCGCCGACCGATGACATCGACATCGGCGATGCCGAACAGTTCGACAGCGGTACGCGCCTGCGCGACACCGACACGCCCGAACGCGAGTTGATGCGCCAGGAACTGGAACAGACGGTGATGAAGGCCGTCAACGCGCTGCCGGAAGAGCTCCGGTCGGCGATCACCCTGCGCGAGGTGGAAGGCCTGAGTACGAGGATATCGCGCAGAAATGGGGTGCCCGATCGGCACCGTGCGTTCACGGATCTTCCGGGCGCGCGAGGCGATCGACACCGAACTCCGGCCGCTGTTGGACATCGGCAGCGCCACCCGTGAGAAGAGCCGCGTATGACCAGTAATCCGTTCAACGAATCGCAGAACCATCAATCGCCCGCCGGGCAGCGCCTGGACCAGCGCCATCGCGAGCAGCTGTCGGCCTGGTCGATGGCGAGCTGGGCCGATGAAGCGCGCTTCCTGCTGCGCCGGATGGAGCACGACCCCGAACTGGCCGGCTGCCAGGAACGCTGGCAGCTGCTGGGCGACGTGATGCGCGGGCAGGCCCCGGCGCTGGCGCCGGCCGGCTTCAGCGCCGCCGTGGCCGCTGCCATCGCCGCCGAGCCGGGCCGGCCGAGCCGCGTCGCCAGGTGCGCCGCAGTGGCTGGCGGGCCTGGGGCGGTGGCGCCGCACTGGCCGCTCGGTGGCGCGTGGCGCTGTTCATGGCGGTGAGAAGCTGAAGGAGGCCACCCGGGCGAGCCGTGGCACCGCAGGTGATCGCCAGCCAGGCGGAACTGGCGCCGGCGCCGACCGACCGGCCCGGTGACCGAGGCGTCGGTGGATACGGCCGCGATGGCCGTGGCGCTGCGCCGGCCGTGGCGATGGCCGCTCCAGCCGCCGCCAGGACACCCGCCGCGCCAGCGCCACCGCACCCAGCAGGCCGCGCGTGCCGCCAGCGCGATGAGGCCCCGCAACGCGCGCCGCAGCGCAGGCGCCGTTGACCCCGACCGTGCCGGCCAATGCCAGCCGCAACATGCCGTTCGGCGAGGTGGGCGGCTGCAGGCGCGACCGTGGCCGCGCTCCAGCCTGGCGCCCGCTGCCGGCGGCGCGCTCAATGCCAGCTTCCCGGCCCATGCCGGGGCGCGGCGTTCTACCCGTTCGAGCCGCGCCTGCAGGACGACCTGCCGGTGCCGCCGCGCCCGCGCGACTGAGCCGGGCGCTTTCACGCCCCGCCGTCGGCGGGGCCCGTATTCCTCTTTCCATCCGATCCGGAGGTTGCCGTCCGATGACTCCCCGACTCCGCACGCAGGCCATGGCCTGCTGGCCCTGACCCTGCCCTGGGGCCTGCGCGCAGGCGCCGGCGCCCGCGCCAGACCCAGCCCGCGCCGCCGCGGCCGCGCGGGCACCGGCCAGCCGCTGGTCACGGCCTGCCCGACTTCACCAACCTGGTCGAACAGGTCGGCCCGGGCGTGGTCAACGTCGACACCACCATCGTCCGCAACAACCGCCAGGCCTCGCGGGCCCGATGGGCGACGATGGCATGCCGGAGTTCTTCCGCCGCTTCTTCGGCCCGGACTTCCCGATGCCGGGGCAGGGCCGGGTGGCCGGACGGTGGCCCCAGCATCAAGGGCCGCGGCATGGGCTCGGGCTTCATCATCTCGCCCGATGGTATGTGCTGACCAACTACCACGTGGTGGCCGATGCCAGCGACGTGAAGGTCAAGCTGGGCGACAGCCGTGAGTTCACCGCCAAGGTGGTGGGCAGCGACCAGCAGTACGACGTCGCGCTGCTGAAGATCGATGGCAAGAACCTGCCGACCGTGCGCGTGGGTGATTCCAACACCCTGAAGCCGGGCCAGTGGGTGTCGCGATCGGCTCGCCGTTCGGCCTCGACCATTCGGTCACCGCGGGCGTGGTCAGCGCCTCGGCCGCAGCACGGGGGGGACCAGCGCTACGTGCCGTTCATCCAGACCGACGTGGCGATCAACCAGGGCAACTCCGGTGGCCCGCTGCTGAACACCCGCGGCGAAGTGGTCGGCATCAACTCGCAGATCTTCTCCGCCTCCGGCGGCTACATGGGCATCAGCTTCGCGATCCCGATCGACCTGGCGATGAGTGCGGTCGAGCAGATCAAGAAGAGTGGCAAGGTCACCCGTGGCCAGCTCGGCGCGGTGGTCGAGCCGATCGACGCTGAAGGCGCAGGGCCTGGGCCTGCCGGACAGCCGTGGCGCCCTGGTCAACCAGATCGTGGCCGGCAGTGCGGCCGCCAAGGCGGGGTGCAGATCGGCGATGTATCCGCTCGGTCAACGGCAGCCCGGTCAACAGCTGGTCCGACCTGCCGCCGCTGATCGGTGCGATGGCACCGGGCAGCCGGTCACCTGGGGTGTCCGCGACGGCAAGCCGCGTGACCTCAGCGCCACCCTGACCGCCTGAGTGAAGACGGCCAGGCCAATGCCCGCGGGCCGGCGGCCGGCGGGGCCGATGCGGCGCCGCAGACCGGCGCCAATGCCCTGCTGGGGCTGGACGTGAGCGACCTGACCGCGCCGCAGCGCAAGCAGTTCGGGCCTGGGCAACGAAGGCGTGCGCATCACCGGGTCAAGGGCCAGGCGCGCGTGATGCCGGGCCTGTCTCCGGGCATGGTGATCCTGCAGGTGGCCCCGGTCGGCAGCGTGGAGCCCTGAACCGGGCGCTGTCCAGCTACAAGAAGGACGACGTGGTGATGCTGCTGGTGCGCACCGGCAACGGCAACAGCGCCTTCGTGGCGGTCAAGGCCGGCCAGTAAGGCCTGCCGGGGCCGCTTCGGCGGCCTGTGTCTCCCCCCGAAGCCCCGCCCGGCGGGGCTTCGGCCGTTCCGGGCGGGCCGATGCGGCCGTTTTCGCGGCGGGGCGCCGGGCATGCGATAATCGACCGTTACCCTGACGACGCCGCGCGCCGCCGCCACCTATGTCTTCTGATTCGATGCGGAACATCCGCAACTTCTCCATCATCGCCCATGTCGACCACGGCAAGTCCACGCTGGCCGACCGCATCATCCAGCTCTGTGGTGGCCTGCAGGCCCGCGAGATGGAAGCGCAGGTGCTCGACTCCAACCCGATCGAGCGCGAGCGTGGCATCACGATCAAGGCGCAGTCGGTGTCGCTGCCGTACCTGGCCAAGGACGGGCAGACCTACCATCTGAACTTCATCGACACCCCCGGCCACGTCGACTTCTCCTATGAAGTCAGCCGCTCGCTGGCCGCCTGCGAAGGCGCGCTGCTGGTGGTCGATGCGGCCAGGGCGTGGAAGCGCAGTCGGTGGCCAACTGCTACACCGCCGTGGAGCAGGGCCTGGAAGTGGTGCCGGTGATCAA";
    //string read = "AAACCCCGTACCGCCACATCTTCAACAACTGCTGCGCGGCTGTACCGCGTCCTGACGGTGCGCACGCCGAGTTTCTGGTAGAGGTTGCGCGATATGCCGTCTTATGGTGTCGCGCCACCGCCAGTTCGCCGCAATCTGTTCGTTGCTGTACCGAATAGATAAGCCTAACACCTGCCCATTCACGCTGGTGAGCGGGCTGGTGCGAAATAAGCTCCGGGAATCCGGATGGTAAGCAGGCGCTCGGACGAAGCTTCGTCAAAATGGGGCGAATTTATGCGTGATGTTGGTTAATCTCTCCGCAGGATGCCGCTGCGCGCGATGCTGCTCCATCTCCGGCAGCGTATTGAGCTGGGTAAGCTTGACGCAATTGTTGCGCCACTGGCTTCTGCCTTCAATTCACAAAATGGGCTGATAAATCCGGTGCGGTTGCGAAGCTGGTAAGGCATCAAGCAATACGCGCTGGGACGTCGTTTTTTTCCGCCCGATGTGCCAAGTAAAGCTGATGGAGTAACAGCAGATTTGCGGTTAAATCGCTCATTAAGCGCAGCCTACGGCATTCTCATTTAACTCTTCAGCACGATTTCCGCAGGCTCGAATTCGCCCAGCAGGATCTGCGCGCGGCGATATTCCGCATTGCCCTTGCAGAAAATGGTTATTGGCAAACCCCGGTTTTGGGCGTATGACGTAACCAGTTTGCGGCCGGTTTTTTGTCGCCGTTAAATACAGCAACAATCACACCGTAAACTTTATCGCGTTAGAGATTCCAGTCGCAGTGATAACGGCCATTGCCTAATAATTTTTCCAGTCGGTTAAGCTGACTACGGGATTATCCCAGGATCGCCCGCGCGCCGCTGGGCACTGCACCAGTAACGTCAGGCACCGCAAACTGTTGTTTGCGGGGAAGGTGACGCACAGCGATATCGCTGCGCCCGGGACCTCCGCTTCGTTCAAGCGCGCCCACGCCCACAGCAATTGCGCGCGGATTACGTACAGAAACTCATGCATCGGCAGTGTTCCAAGATGCTGCTCCTTTTGATCAATGGAAGGCTCTCTCCTGCGTTTTCCCATGCGGCCTGCAGAATCGCCCTGGCGAACTGAATTTCCGCTTGTTGAATTAAGCTCCA";
    //string ref = "AAACCCCGTACCCCATCATCTTCAGCAACTGCTGCGCGTGCTGTACCGCGTCCTGACGGTGCGCGACGCCGAGTTTCTGGTAGAGGTTGCGGATATGCGTCTTAATGGTGGTCGCGGCCACCGCCAGTTCGCCGGCAATCTGTTCGTTGCTGTAACCGGAATAGATAAGCCCTAACACCTGCCATTCACGCTGGGTGAGCGGGCTGGTGCGAATAAGCTCCGGGACATCCGGATGGTTAAGCAGGCGCTCGACGAAGCCTTCGTCAAAATGGGCGAATTTATGCCGGTGATGTTGGTTAATCTCTCGCAGGATGCGCTGCGCGCGATGCTGCTCCATCTCCGGCAGCGTATTGAGCTGGATAAGCTGACGCAATTGTTGCGCCATGGCTTCGCCTTCAATCACAAAATGGCTGATAAATCCGGTGCGGTTGGCGAGCTGTAAGGCATCAAGCAATACGCGCTGGGCGTCGTTTTTCCGCCCGGATTGCCAGTAAAGCTGATTGAGTAACAGCAGATTGCGGTTTAAATCGCTCATTAAGCGCAGGCTACGGGCATTCTCATTTAACTCTTCCAGCACGATTTCCGCAGGCTCGAATTCGCCCAGCAGGATCTGCGCGCGGGCGATATTCCGCCATTGCCCTTGCAGAAAATGGTTATTGGCAAACGCCGGTTTGGGCGTATGACGTAACCAGTTGGCGGCGGATTTTTTGTCGCCGGTTAACTGCCAGTAAATGACACGCACCTTATCGGCGTTAGAGATCCAGTCGCAGTGATAACGGCCATTGCCTAATAAATTTTCCAGTCGGTTAAGCTGACTACGGGCATTATCCAGATCGCCGCGCGCCAGCGAGCACTGCACCAGTAACGTCAGGCACTGCAACTGTTGTTGCGGCTGGAAGGTGGACAGCACAGCGATACCGCTGCGCGCCGAGGCCTCCGCTTCGTCAAGGCGCGCCCACGCCCACAGCAATTGCGCGCGGATACGTACCAGAAACTCATGCATCGGCAGTTGTTCCAGATGTTGCTCTTTGATCAACTGGAAGGCTCTCTCCTGCGTTTCCCATGCGGCCTGCAGGAATCCCTGGGCGAACTGAATTTCGCTTTGTTGAATTAAGCTCCA";
    
    StripedSmithWaterman::Alignment result;
    Alignment aligner;
    aligner.local_align(read, ref, result);
    
    cout << result.cigar_string << endl;
    AlignCoderSNV aligncodersnv;
    vector<int> encode_data;
    aligncodersnv.encode(result, read, ref, 0, encode_data);
    cout << encode_data << endl;
    
}

TEST_CASE("test AlignCoderSNV::encode (ssw) for a whole .m5 file")
{
    string m5_file = "../data/MSSA_61_forward.m5";
    string encode_file = "../results/MSSA_61_forward.encode.ssw";
    ofstream fs_encode; open_outfile(fs_encode, encode_file);
    
    AlignReaderM5 alignreader;
    alignreader.open(m5_file);
    Align align;
    int nline = 0;
    while(alignreader.readline(align)){
        ++nline;
        
        // expections
        int alen = (int) align.matchPattern.size();
        if ( !(align.qAlignedSeq.size()==alen && align.tAlignedSeq.size()==alen) )
            throw runtime_error("incorrect match patter in line " + to_string(nline));
        if (align.qStrand != '+')
            throw runtime_error("qStrand should be + in line " + to_string(nline));
        
        // remove indel from alignment pattern in m5 file
        string qSeq = rm_indel_from_seq(align.qAlignedSeq);
        string tSeq = rm_indel_from_seq(align.tAlignedSeq);
        
        // align qSeq to tSeq by ssw
        StripedSmithWaterman::Alignment result;
        Alignment aligner;
        aligner.local_align(qSeq, tSeq, result);

        // encode
        AlignCoderSNV aligncodersnv;
        vector<int> encode_data;
        aligncodersnv.encode(result, qSeq, tSeq, align.tStart, encode_data);

        // print encode
        for (int i=0; i<(int)encode_data.size(); ++i)
            fs_encode << encode_data[i] << "\t";
        fs_encode << endl;
    }
    
    alignreader.close();
    fs_encode.close();
    
}




