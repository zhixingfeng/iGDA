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
    
    // check if alignment is solf clipped
    for (auto i = 0; i < seqan::length(record.cigar); ++i){
        if (record.cigar[i].operation == 'S' || record.cigar[i].operation == 's')
            throw runtime_error("the sam file is soft clipped, try hard clip");
    }
    
    if (seqan::length(this->ref_ids) == 0 || seqan::length(this->ref_seqs) == 0)
        return true;
    // get alignment if reference file is provided
    typedef seqan::Align<seqan::Dna5String> TAlign;
    typedef seqan::Row<TAlign>::Type TRow;
    typedef seqan::Iterator<TRow>::Type TRowIterator;
    
    TAlign cur_align;
    
    if (record.rID != seqan::BamAlignmentRecord::INVALID_REFID)
        bamRecordToAlignment(cur_align, ref_seqs[record.rID], record);

    //cout << cur_align << endl;
    TRow row_ref = seqan::row(cur_align,0);
    TRow row_read = seqan::row(cur_align,1);
    if (seqan::length(row_ref) != seqan::length(row_read)){
        cout << "row_ref = " << row_ref << endl;
        cout << "row_read = " << row_read << endl;
        cout << "seqan::length(row_ref) = " << seqan::length(row_ref) << endl;
        cout << "seqan::length(row_read) = " << seqan::length(row_read) << endl;
        throw runtime_error("AlignReaderSam::readline: seqan::length(row_ref) != seqan::length(row_read)");
    }
    align.tAlignedSeq.resize(seqan::length(row_ref));
    align.qAlignedSeq.resize(seqan::length(row_read));
    align.matchPattern.resize(seqan::length(row_ref));
    
    TRowIterator it_row_ref = begin(row_ref);
    TRowIterator it_end = end(row_ref);
    TRowIterator it_row_read = begin(row_read);
    
    int i = 0;
    for (; it_row_ref != it_end; ++it_row_ref, ++it_row_read){
        if (!seqan::isGap(it_row_ref)){
            align.tAlignedSeq[i] = seqan::Dna5(*it_row_ref);
            if (!seqan::isGap(it_row_read)){
                align.qAlignedSeq[i] = seqan::Dna5(*it_row_read);
                if (align.tAlignedSeq[i] == align.qAlignedSeq[i])
                    align.matchPattern[i] = '|';
                else
                    align.matchPattern[i] = '*';
            }else{
                align.qAlignedSeq[i] = '-';
                align.matchPattern[i] = '*';
            }
        }else{
            align.tAlignedSeq[i] = '-';
            if (!seqan::isGap(it_row_read)){
                align.qAlignedSeq[i] = seqan::Dna5(*it_row_read);
                align.matchPattern[i] = '*';
            }else{
                throw runtime_error("AlignReaderSam::readline: gap vs gap!");
            }
        }
        ++i;
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








