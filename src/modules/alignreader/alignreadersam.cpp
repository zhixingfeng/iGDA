//
//  alignreadersam.cpp
//  iGDA
//
//  Created by Zhixing Feng on 6/29/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "alignreadersam.h"
#include "../../misc/basic.h"

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
    
    if(seqan::length(record.qual) != seqan::length(record.seq) && seqan::length(record.qual) > 0)
        throw runtime_error("align.qv.size() != seqan::length(record.seq)");
    //cout << seqan::length(record.qual) << endl;
    int64_t seq_shift = 0;
    int64_t ref_shift = 0;
    align.qAlignedSeq.clear();
    align.tAlignedSeq.clear();
    align.matchPattern.clear();
    align.qv.clear();
    align.qv_locus.clear();
    for (auto i = 0; i < seqan::length(record.cigar); ++i){
        switch (record.cigar[i].operation) {
            case 'M':
                // match or mismatch
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( record.seq[seq_shift] );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    
                    if (seqan::length(record.qual) > 0){
                        align.qv.push_back((uint8_t)record.qual[seq_shift] - uint8_t(33));
                        align.qv_locus.push_back(record.beginPos + ref_shift);
                    }
                    
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
                    
                    if (seqan::length(record.qual) > 0){
                        align.qv.push_back((uint8_t)record.qual[seq_shift] - uint8_t(33));
                        align.qv_locus.push_back(record.beginPos + ref_shift);
                    }
                    
                    align.matchPattern.push_back('*');
                    ++seq_shift;
                }
                break;
            case 'D':
                // deletion
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( '-' );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    
                    if (seqan::length(record.qual) > 0){
                        align.qv.push_back(0);
                        align.qv_locus.push_back(record.beginPos + ref_shift);
                    }
                    
                    align.matchPattern.push_back('*');
                    ++ref_shift;
                }
                break;
            case 'N':
                // skip, same as D
                for (auto j = 0; j < record.cigar[i].count; ++j){
                    align.qAlignedSeq.push_back( '-' );
                    align.tAlignedSeq.push_back( cur_refseq[record.beginPos + ref_shift] );
                    
                    if (seqan::length(record.qual) > 0){
                        align.qv.push_back(0);
                        align.qv_locus.push_back(record.beginPos + ref_shift);
                    }
                    
                    align.matchPattern.push_back('*');
                    ++ref_shift;
                }
                break;
            case 'S':
                // soft clip, ignore
                for (auto j = 0; j < record.cigar[i].count; ++j)
                    ++seq_shift;
                if (i == 0)
                    align.qStart = record.cigar[i].count;
                break;
            case 'H':
                // hard clip, ignore
                if (i == 0)
                    align.qStart = record.cigar[i].count;
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
                    
                    if (seqan::length(record.qual) > 0){
                        align.qv.push_back((uint8_t)record.qual[seq_shift] - uint8_t(33));
                        align.qv_locus.push_back(record.beginPos + ref_shift);
                    }
                    
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
                    
                    if (seqan::length(record.qual) > 0){
                        align.qv.push_back((uint8_t)record.qual[seq_shift] - uint8_t(33));
                        align.qv_locus.push_back(record.beginPos + ref_shift);
                    }
                    
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
   
    if (align.qAlignedSeq.size() != align.tAlignedSeq.size())
        throw runtime_error("align.qAlignedSeq.size() != align.tAlignedSeq.size()");
    if ( !(align.qAlignedSeq.size() ==  align.qv.size() && align.qAlignedSeq.size() ==  align.qv_locus.size()) &&
        seqan::length(record.qual) > 0)
        throw runtime_error("incompatible alignment pattern with quality score");
    //cout << align.qAlignedSeq.size() << endl;
    return true;
}

bool AlignReaderSam::close() {
    seqan::close(bamFileIn);
    return true;
}


bool AlignReaderSam::read(string filename, vector<Align> &align_vec)
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
    
    for (auto i = 0; i < seqan::length(ref_ids); ++i){
        stringstream ss;
        ss << ref_ids[i];
        vector<string> buf = split(ss.str(), ' ');
        if (buf.size() == 0)
            throw runtime_error("AlignReaderSam::getref(): buf.size() == 0");
        ref_ids[i] = buf[0];
        //cout << ref_ids[i] << endl;
    }
    seqan::close(seqFileIn);
    return true;
}

bool AlignReaderSam::samtom5(string sam_file, string ref_file, string m5_file, int min_len)
{
    this->getref(ref_file);
    this->open(sam_file);
    ofstream fs_m5_file;
    open_outfile(fs_m5_file, m5_file);
    Align align;
    while(this->readline(align)){
        if (align.tEnd - align.tStart + 1 < min_len)
            continue;
        fs_m5_file << align.qName << ' ' << '0' << ' ' << align.qStart << ' ' << '0' << ' ' << align.qStrand << ' ' << align.tName << ' ' << '0' << ' ';
        fs_m5_file << align.tStart << ' ' << align.tEnd + 1 << ' ' << align.tStrand << ' ';
        fs_m5_file << '0' << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ';
        fs_m5_file << align.mapQV << ' ' << align.qAlignedSeq << ' ' << align.matchPattern << ' ' << align.tAlignedSeq << endl;
    }
    this->close();
    fs_m5_file.close();
    return true;
}

bool AlignReaderSam::samtom5qv(string sam_file, string ref_file, string m5qv_file, int min_len)
{
    this->getref(ref_file);
    this->open(sam_file);
    ofstream fs_m5qv_file;
    open_outfile(fs_m5qv_file, m5qv_file);
    Align align;
    while(this->readline(align)){
        if (align.tEnd - align.tStart + 1 < min_len)
            continue;
        fs_m5qv_file << align.qName << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ' << align.qStrand << ' ' << align.tName << ' ' << '0' << ' ';
        fs_m5qv_file << align.tStart << ' ' << align.tEnd + 1 << ' ' << align.tStrand << ' ';
        fs_m5qv_file << '0' << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ' << '0' << ' ';
        fs_m5qv_file << align.mapQV << ' ' << align.qAlignedSeq << ' ' << align.matchPattern << ' ' << align.tAlignedSeq << ' ';
        for (auto i = 0; i < align.qv.size(); ++i){
            if (i < align.qv.size() - 1){
                fs_m5qv_file << uint32_t(align.qv[i]) << ',';
            }else{
                fs_m5qv_file << uint32_t(align.qv[i]);
            }
        }
        fs_m5qv_file << ' ';
        fs_m5qv_file << align.qv_locus << endl;
    }
    this->close();
    fs_m5qv_file.close();
    return true;
}

bool AlignReaderSam::getchrrange(string sam_file, string ref_file, string out_file)
{
    this->getref(ref_file);
    this->open(sam_file);
    Align align;
    unordered_map<string, pair<size_t, size_t> > chrrange;
    size_t n_read = 0;
    while(this->readline(align)){
        auto it = chrrange.find(align.tName);
        if (it == chrrange.end()){
            chrrange[align.tName].first = align.tStart;
            chrrange[align.tName].second = align.tEnd;
        }else{
            if (align.tStart < it->second.first) it->second.first = align.tStart;
            if (align.tEnd > it->second.second) it->second.second = align.tEnd;
        }
        ++n_read;
        if (n_read % 1000 == 0) cout << n_read << endl;
    }
    this->close();
    
    ofstream fs_out_file;
    open_outfile(fs_out_file, out_file);
    for (auto it = chrrange.begin(); it != chrrange.end(); ++it){
        fs_out_file << it->first << '\t' << it->second.first << '\t' << it->second.second << endl;
        //cout << it->first << '\t' << it->second.first << '\t' << it->second.second << endl;
    }
    fs_out_file.close();
    
    return true;
}


