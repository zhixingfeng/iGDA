//
//  test_AlignCoder.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/5/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
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



TEST_CASE("test AlignCoderSNV::encode_ssw")
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
    aligncodersnv.encode_ssw(result, read, ref, 0, encode_data);
    cout << encode_data << endl;
    
}


