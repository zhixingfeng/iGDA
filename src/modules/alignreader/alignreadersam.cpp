//
//  alignreadersam.cpp
//  iGDA
//
//  Created by Zhixing Feng on 6/29/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "alignreadersam.h"


bool AlignReaderSam::open(string filename) {
    // Open input file, BamFileIn can read SAM and BAM files.
    if (!seqan::open(bamFileIn, filename.c_str()))
    {
        throw runtime_error("AlignReaderSam::open: fail to open " + filename);
        return false;
    }

    // check reference in header consistent with loaded fasta file
    seqan::  BamHeader header;
    readHeader(header, bamFileIn);
    
    typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;
    TBamContext const & bamContext = context(bamFileIn);
    
    if (seqan::length(ref_ids) == 0 || seqan::length(ref_seqs) == 0)
        return true;
    
    // check number of contigs
    if (seqan::length(seqan::contigNames(bamContext)) != seqan::length(ref_ids) || seqan::length(seqan::contigNames(bamContext)) != seqan::length(ref_seqs))
        throw runtime_error("ref in sam is not consistent with loaded fasta file");
    
    // check contig names and length of sequences
    for (unsigned i = 0; i < seqan::length(seqan::contigNames(bamContext)); ++i){
        if (seqan::contigNames(bamContext)[i] != ref_ids[i])
            throw runtime_error("ref name in sam and loaded fasta do not match");
        if (seqan::contigLengths(bamContext)[i] != seqan::length(ref_seqs[i]))
            throw runtime_error("ref sequence in sam and loaded fasta do not match");
    }

    return true;
}

bool AlignReaderSam::readline(Align &align) {
    if (atEnd(bamFileIn))
        return false;
    seqan::BamAlignmentRecord record;
    readRecord(record, bamFileIn);
    
    align.qName = toCString(record.qName);
    align.qStrand = '+';
    align.tStart = record.beginPos;
    align.tEnd = align.tStart + getAlignmentLengthInRef(record) - 1;
    align.tStrand = hasFlagRC(record) ? '-' : '+';
    align.mapQV = record.mapQ;
    
    // get alignment of the current record
    if (record.rID == seqan::BamAlignmentRecord::INVALID_REFID)
        throw runtime_error("AlignReaderSam::readline: invalid refid");
    
    if (record.rID >= seqan::length(this->ref_ids)){
        throw runtime_error("AlignReaderSam::readline: unavailable refID " + to_string(record.rID) );
    }
    
    seqan::Dna5String cur_refseq = ref_seqs[record.rID];
    align.tName = seqan::toCString(ref_ids[record.rID]);
    int64_t seq_shift = 0;
    int64_t ref_shift = 0;
    align.qAlignedSeq.clear();
    align.tAlignedSeq.clear();
    align.matchPattern.clear();
    for (auto i = 0; i < seqan::length(record.cigar); ++i){
        switch (record.cigar[i].operation) {
            case 'M':
                // match or mismatch
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( record.seq[seq_shift] );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    if (record.seq[seq_shift] == cur_refseq[record.beginPos + ref_shift]){
                        align.matchPattern.push_back('|');
                    }else{
                        align.matchPattern.push_back('*');
                    }
                    ++seq_shift;
                    ++ref_shift;
                }
                break;
                
            case 'I':
                // insertion
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( record.seq[seq_shift] );
                    align.tAlignedSeq.push_back( '-' );
                    align.matchPattern.push_back('*');
                    ++seq_shift;
                }
                break;
            case 'D':
                // deletion
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( '-' );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    align.matchPattern.push_back('*');
                    ++ref_shift;
                }
                break;
            case 'N':
                // skip, same as D
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( '-' );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    align.matchPattern.push_back('*');
                    ++ref_shift;
                }
                break;
            case 'S':
                // soft clip, ignore
                for (auto j = 0; j < record.cigar[i].count; ++j)
                    ++seq_shift;
                break;
            case 'H':
                // hard clip, ignore and do nothing
                break;
            case 'P':
                // forbid padding sam
                throw runtime_error("AlignReaderSam::readline: padding SAM is not allowed");
                break;
            case '=':
                // match
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( record.seq[seq_shift] );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    if (record.seq[seq_shift] == cur_refseq[record.beginPos + ref_shift]){
                        align.matchPattern.push_back('|');
                    }else{
                        throw runtime_error("AlignReaderSam::readline: cigar is = but seq and ref are different");
                    }
                    ++seq_shift;
                    ++ref_shift;
                }
                break;
            case 'X':
                // mismatch
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( record.seq[seq_shift] );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    if (record.seq[seq_shift] == cur_refseq[record.beginPos + ref_shift]){
                        throw runtime_error("AlignReaderSam::readline: cigar is = but seq and ref are different");
                    }else{
                        align.matchPattern.push_back('*');
                    }
                    ++seq_shift;
                    ++ref_shift;
                }
                break;
                
            default:
                throw runtime_error("AlignReaderSam::readline: unknown cigar type");
                break;
        }
        
    }
   
    
    return true;
}

bool AlignReaderSam::close() {

    return true;
}


bool AlignReaderSam::read(string filename, stxxl::vector<Align> &align_vec)
{
        
    return true;
}


 bool AlignReaderSam::getref(string filename)
{
    seqan::SeqFileIn seqFileIn;
    if (!seqan::open(seqFileIn, filename.c_str()))
    {
        throw runtime_error("AlignReaderSam::open: fail to open " + filename);
        return false;
    }
    
    try{
        seqan::readRecords(ref_ids, ref_seqs, seqFileIn);
    }catch (seqan::Exception const & e){
        std::cout << "ERROR: " << e.what() << std::endl;
        return false;
    }
    
    return true;
}

bool AlignReaderSam::samtom5(string sam_file, string ref_file, string m5_file)
{
    this->getref(ref_file);
    this->open(sam_file);
    ofstream fs_m5_file;
    open_outfile(fs_m5_file, m5_file);
    Align align;
    while(this->readline(align)){
        fs_m5_file << align.qName << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ' << align.qStrand << ' ' << align.tName << ' ' << '0' << ' ';
        fs_m5_file << align.tStart << ' ' << align.tEnd + 1 << ' ' << align.tStrand << ' ';
        fs_m5_file << '0' << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ';
        fs_m5_file << align.mapQV << ' ' << align.qAlignedSeq << ' ' << align.matchPattern << ' ' << align.tAlignedSeq << endl;
    }
    this->close();
    fs_m5_file.close();
    return true;
}






