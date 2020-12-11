//
//  var2vcf.cpp
//  iGDA
//
//  Created by Zhixing Feng on 12/10/20.
//  Copyright © 2020 Zhixing Feng. All rights reserved.
//

#include "var2vcf.h"
#include "../modules/alignreader/alignreadersam.h"
#include "../../include/utils.h"


void var2vcf(string var_file, string ref_file, string vcf_file)
{
    AlignReaderSam sam_reader;
    sam_reader.getref(ref_file);
    if (seqan::length(sam_reader.ref_ids) != 1){
        throw runtime_error("var2vcf(): seqan::length(sam_reader.ref_ids) != 1");
    }
    
    const seqan::Dna5String &cur_refseq = sam_reader.ref_seqs[0];
    const string cur_refid = seqan::toCString(sam_reader.ref_ids[0]);
    std::time_t cur_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ofstream fs_vcffile;
    open_outfile(fs_vcffile, vcf_file);
    // print vcf information
    fs_vcffile << "##fileformat=VCFv4.2" << endl;
    fs_vcffile << "##fileDate=" << std::ctime(&cur_time);
    fs_vcffile << "##source=igda" << endl;
    fs_vcffile << "##reference=" << ref_file << endl;
    fs_vcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << endl;
    
    // print SNVs
    ifstream fs_varfile;
    open_infile(fs_varfile, var_file);
    int64_t pre_pos = -1;
    vector<string> buf_vec;
    size_t nr = 1;
    string pre_altbase = "";
    char pre_refbase = '@';
    while (true) {
        // read line of var file
        string buf;
        getline(fs_varfile, buf);
        if (fs_varfile.eof())
            break;
        buf_vec = split(buf, '\t');
        int64_t cur_pos = stoll(buf_vec[0].c_str());
        if (buf_vec[1].size() != 1){
            throw runtime_error("var2vcf(): buf_vec[1].size() != 1");
        }
        
        // check the data
        string altbase = buf_vec[1];
        char refbase = cur_refseq[cur_pos];
        if (altbase[0] == refbase){
            throw runtime_error("var2vcf(): altbase == refbase");
        }
        
        // assuming the var file is sorted!!
        if (cur_pos == pre_pos){
            pre_altbase = pre_altbase + "," + altbase;
        }else{
            if (nr > 1){
                fs_vcffile << cur_refid << '\t' << pre_pos + 1 << "\t.\t" << pre_refbase << '\t' << pre_altbase << "\t100\tPASS\t." << endl;
            }
            pre_pos = cur_pos;
            pre_refbase = refbase;
            pre_altbase = altbase;
        }
        ++nr;
    }
    fs_vcffile << cur_refid << '\t' << pre_pos + 1 << "\t.\t" << pre_refbase << '\t' << pre_altbase << "\t100\tPASS\t." << endl;
    fs_varfile.close();
    fs_vcffile.close();
}


