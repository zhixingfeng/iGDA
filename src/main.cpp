//
//  main.cpp
//  iGDA
//
//  Created by Zhixing Feng on 8/4/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//
#define CATCH_CONFIG_RUNNER
#include "../include/catch.hpp"

#include <headers.h>
#include "../tools/tools.h"
#include "../src/modules/modules.h"
#include "../src/modules/dforest/dforestsnvmax.h"
#include "../src/modules/dforest/dforestsnvstxxl.h"
#include "../src/modules/errormodel/errormodelsnv.h"
#include "../src/modules/assemble/assembler.h"
#include "../src/modules/detectsingle/detectsinglesnv.h"

#include "./misc/misc.h"

#ifdef _UNITTEST

using namespace boost::filesystem;

int main(int argc, char* argv[]) {
    if (!exists("../results"))
        create_directory("../results");
    int result = Catch::Session().run( argc, argv );
    return result;
}

#else

using namespace TCLAP;

void print_usage()
{
    cout << "igda [command]" << endl;
    cout << "command = encode : binary coding SNVs of each aligned read." << endl;
    cout << "          cmpreads : pairwise comparison of encoded reads."<< endl;
    cout << "          dforest : detecting SNVs by dforest algorithm (unable to detect SNVs if no reads cover multiple real SNVs)." << endl;
    cout << "          contexteffect : pileup reads and get context effect." << endl;
    cout << "          detectsingle : detecting SNVs locus by locus (lower sensitivity compared to dforest but able to detect SNVs if no reads cover multiple real SNVs)." << endl;
    cout << "          rdim : removing undetected SNVs from encoded reads." << endl;
    cout << "          ann : clustering reads by adaptive nearest neighbor clustering." << endl;
    
    //cout << "command = samtofa, bamtofa, m5tofa, encode, cmpreads, bin2txt, txt2bin, dforest, sort, filter, contexteffect, merge, mergeall, dist, pileup_var, pileup_reads, samtom5" << endl;

}

int main(int argc, const char * argv[])
{
    try{
        if (argc == 1) { print_usage(); return 0; }
        
        // load arguments
        string cmdname = "igda "; cmdname += argv[1];
        vector<string> argv2;
        argv2.push_back(cmdname);
        for (int i=2; i<argc; i++){
            argv2.push_back(argv[i]);
        }
        
        CmdLine cmd("iGDA", ' ', "0.1");
        
        // sam to fasta
        if (strcmp(argv[1], "samtofa")==0) {
            UnlabeledValueArg<string> samfileArg("samfile", "path of sam file", true, "", "samfile", cmd);
            UnlabeledValueArg<string> fafileArg("fafile", "path of fa file", true, "", "fafile", cmd);
            
            cmd.parse(argv2);
            
            string shell_cmd = "samtools view " + samfileArg.getValue() +
            " |  awk \'{print \">\"$1; print $10}\' > " +  fafileArg.getValue();
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            return 0;
        }

        
        // bam to fasta
        if (strcmp(argv[1], "bamtofa")==0) {
            UnlabeledValueArg<string> bamfileArg("bamfile", "path of bam file", true, "", "bamfile", cmd);
            UnlabeledValueArg<string> fafileArg("fafile", "path of fa file", true, "", "fafile", cmd);
            UnlabeledValueArg<string> chrArg("chr", "chromosome", false, "", "chr", cmd);
            
            cmd.parse(argv2);
            
            string shell_cmd = "samtools view " + bamfileArg.getValue() + " " + chrArg.getValue() + 
                                " | awk \'{print \">\"$1; print $10}\' > " +  fafileArg.getValue();
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            return 0;
        }

        
        // m5 to fasta
        if (strcmp(argv[1], "m5tofa")==0) {
            UnlabeledValueArg<string> m5fileArg("m5file", "path of m5 file", true, "", "m5file", cmd);
            UnlabeledValueArg<string> fafileArg("fafile", "path of fa file", true, "", "fafile", cmd);
            cmd.parse(argv2);
            
            m5tofa(m5fileArg.getValue(), fafileArg.getValue());
            return 0;
        }
        
        // encode alignment file
        if (strcmp(argv[1], "encode")==0){
            ValueArg<int> mArg("m","method","method for encoding. 0: full, 1: SNV", false , 1, "method", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of alignment file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", false, "", "reffile", cmd);
            SwitchArg ism5Arg("t", "m5", "is use m5 format", cmd, false);
            
            cmd.parse(argv2);
            
            // set alignreader
            AlignReaderM5 alignreaderm5;
            AlignReaderSam alignreadersam;
            
            // select encoder
            AlignCoderSNV aligncodersnv;
            AlignCoder *p_aligncoder = NULL;
            
            switch (mArg.getValue()){
                case 0:
                    return 0;
                    break;
                case 1:
                    p_aligncoder = &aligncodersnv;
                    break;
                default:
                    cerr << "Error: invalid argument. -m should be 0 or 1" << endl;
                    return 1;
            }
            
            // setup reader
            if (ism5Arg.getValue()){
                p_aligncoder->setAlignReader(&alignreaderm5);
            }else{
                if (reffileArg.getValue() == "")
                    throw runtime_error("use sam file as input but not reference file provided.");
                alignreadersam.getref(reffileArg.getValue());
                p_aligncoder->setAlignReader(&alignreadersam);
            }
            p_aligncoder->encode(alignfileArg.getValue(), outfileArg.getValue());
            
            return 0;
        }
        
        // pairwise compare reads
        if (strcmp(argv[1], "cmpreads")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            ValueArg<int> topnArg("p","topn","select top n candidates of each reads, default: 20", false , 20, "topn", cmd);
            ValueArg<double> overlapArg("l","overlap","minimal overlap of reads, default: 0.5", false , 0.5, "overlap", cmd);
            SwitchArg istextArg("t", "text", "is output text file", cmd, false);
            SwitchArg isreadIDArg("r", "readID", "is print reads ID", cmd, false);
            SwitchArg isnocondprobArg("c", "condprob", "not use conditional probability", cmd, false);
            SwitchArg isdiffArg("d", "diff", "is calculate difference between reads", cmd, false);
            
            SwitchArg islegacyArg("y", "legacy", "is use the legacy algorithm", cmd, false);
            
            cmd.parse(argv2);
            
            if (isdiffArg.getValue()){
                cmpreads_topn_diff(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                                   topnArg.getValue(), overlapArg.getValue());
            }else{
                if (!islegacyArg.getValue()){
                    cmpreads_topn(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                              topnArg.getValue(), overlapArg.getValue(), true, !istextArg.getValue(),
                              isreadIDArg.getValue(), !isnocondprobArg.getValue());
                }else{
                    cmpreads_topn_legacy(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                                  topnArg.getValue(), overlapArg.getValue(), true, !istextArg.getValue(),
                                  isreadIDArg.getValue(), !isnocondprobArg.getValue());

                }
                
            }
            
            
        }
        
        // convert binary cmpreadsfile to text
        if (strcmp(argv[1], "bin2txt") == 0) {
            UnlabeledValueArg<string> binfileArg("binfile", "binary cmpreads file", true, "", "binfile", cmd);
            UnlabeledValueArg<string> txtfileArg("txtfile", "text cmpreads file", true, "", "txtfile", cmd);
            SwitchArg isreadIDArg("r", "readID", "is use reads ID", cmd, false);
            SwitchArg isdiffArg("d", "diff", "is use differences", cmd, false);
            
            cmd.parse(argv2);
            
            cmpreads_bin2txt(binfileArg.getValue(), txtfileArg.getValue(), isreadIDArg.getValue(), isdiffArg.getValue());
        }
        
        // convert text cmpreadsfile to binary
        if (strcmp(argv[1], "txt2bin")==0) {
            UnlabeledValueArg<string> txtfileArg("txtfile", "text cmpreads file", true, "", "txtfile", cmd);
            UnlabeledValueArg<string> binfileArg("binfile", "binary cmpreads file", true, "", "binfile", cmd);
            SwitchArg isreadIDArg("r", "readID", "is use reads ID", cmd, false);
            
            cmd.parse(argv2);
            
            cmpreads_txt2bin(txtfileArg.getValue(), binfileArg.getValue(), isreadIDArg.getValue());
        }
        
        // dforest algorithm
        if (strcmp(argv[1], "dforest")==0){
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> cmpreadsfileArg("cmpreadsfile", "path of cmpreads file", true, "", "cmpreadsfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            UnlabeledValueArg<string> tmpdirArg("tmpdir", "temporary directory", true, "", "tmpdir", cmd);
            
            ValueArg<int> minreadsArg("r","minreads","minimal number of reads in a node, default: 8", false , 8, "minreads", cmd);
            ValueArg<int> maxdepthArg("d","maxdepth","maximal depth of a tree, default: 5", false , 5, "maxdepth", cmd);
            ValueArg<double> minfreqArg("f","minfreq","minimal frequency: 0.0", false , 0.0, "minfreq", cmd);
            ValueArg<int> nthreadArg("n","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            
            SwitchArg islegacyArg("l", "legacy", "use the legacy algorithm (no stxxl) to run dforest", cmd, false);
            SwitchArg isinterArg("i", "intermediate", "output intermediate results", cmd, false);
            
            cmd.parse(argv2);
            cout << "minreads = " << minreadsArg.getValue() << endl;
            cout << "maxdepth = " << maxdepthArg.getValue() << endl;
            cout << "minfreq = " << minfreqArg.getValue() << endl;
            cout << "nthread = " << nthreadArg.getValue() << endl;
            
            AlignReaderM5 alignreaderm5;
            AlignReaderSam alignreadersam;
            AlignCoderSNV aligncoder;
            
            AlignReader *p_alignreader = NULL;
            string alignfile = alignfileArg.getValue();
            if (alignfile.substr(alignfile.size()-4, 4) == ".sam"){
                cout << "use sam file " << endl;
                p_alignreader = &alignreadersam;
            }else{
                p_alignreader = &alignreaderm5;
            }
            DForestSNVMax forestsnvmax(p_alignreader, &aligncoder);
            DForestSNVSTXXL forestsnvstxxl(p_alignreader, &aligncoder);
            
            DForest *ptr_forest;
            if (islegacyArg.getValue())
                ptr_forest = &forestsnvmax;
            else
                ptr_forest = &forestsnvstxxl;
            
            string shell_cmd = "mkdir -p " + tmpdirArg.getValue();
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            
            ptr_forest->run(encodefileArg.getValue(), alignfileArg.getValue(), cmpreadsfileArg.getValue(), outfileArg.getValue(), tmpdirArg.getValue(),
                            minreadsArg.getValue(), maxdepthArg.getValue(), nthreadArg.getValue(), minfreqArg.getValue(), isinterArg.getValue());
        }
        
        // sort output
        if (strcmp(argv[1], "sort")==0) {
            UnlabeledValueArg<string> fileArg("outputfile", "", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            string sortfile = fileArg.getValue() + ".sorted";
            string maxfile = fileArg.getValue() + ".max";
            
            string shell_cmd;
            
            shell_cmd = "sort -s -k1,1n -k3,3nr " + fileArg.getValue() + " > " + sortfile;
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            
            shell_cmd = "sort -s -u -k1,1n " + sortfile + " > " + maxfile;
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
        }
        
        // get all variants
        if (strcmp(argv[1], "getvar")==0) {
            UnlabeledValueArg<string> dforestfileArg("dforestfile", "path of dforest file", true, "", "dforestfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            ValueArg<double> mincondprobArg("c","condprob","minimal conditional probability, default: 0.7", false , 0.7, "condprob", cmd);
            
            cmd.parse(argv2);
            
            Assembler assembler;
            assembler.get_variants(dforestfileArg.getValue(), outfileArg.getValue(), mincondprobArg.getValue());
        }
        
        // remove variants out of var_file to reduce dimension
        if (strcmp(argv[1], "rdim")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> varfileArg("varfile", "path of var file", true, "", "varfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of out file", true, "", "outfile", cmd);
            cmd.parse(argv2);
            
            Assembler assembler;
            assembler.reduce_dim(encodefileArg.getValue(), varfileArg.getValue(), outfileArg.getValue());
        }
        // filter output
        if (strcmp(argv[1], "filter")==0) {
            UnlabeledValueArg<string> dforestfileArg("dforestfile", "path of dforest output file", true, "", "dforestfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            ValueArg<double> minfreqArg("f","minfreq","minimal frequency: 0.0", false , 0.0, "minfreq", cmd);
            
            cmd.parse(argv2);
            
            string shell_cmd = "awk -F'\\t' '{if ($3>=" + to_string(minfreqArg.getValue()) + ") print}' " + 
                                dforestfileArg.getValue() + " > " + outfileArg.getValue();
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            
        }
        
        // learn context effect
        if (strcmp(argv[1], "contexteffect")==0) {
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of alignment file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outprefixArg("outprefix", "prefix of output files", true, "", "outprefix", cmd);
            
            ValueArg<int> leftlenArg("l","leftlen","length of the upstream context, default: 1", false , 1, "leftlen", cmd);
            ValueArg<int> rightlenArg("r","rightlen","length of the downstream context, default: 1", false , 1, "rightlen", cmd);
            
            cmd.parse(argv2);
            ErrorModelSNV errormodel;
            errormodel.set_context_size(leftlenArg.getValue(), rightlenArg.getValue());
            errormodel.learn(alignfileArg.getValue(), outprefixArg.getValue());
        }
        
        // merge context effect
        if (strcmp(argv[1], "merge")==0) {
            cmd.parse(argv2);
            vector<string> context_files;
            string buf;
            while (getline(cin, buf)) {
                context_files.push_back(buf);
            }
                        
            ErrorModelSNV errormodel;
            errormodel.merge(context_files);
        }
        
        // merge context effect all
        if (strcmp(argv[1], "mergeall")==0) {
            cmd.parse(argv2);
            vector<string> context_files;
            string buf;
            while (getline(cin, buf)) {
                context_files.push_back(buf);
            }
            
            ErrorModelSNV errormodel;
            errormodel.merge_all(context_files);
        }

        
        // calculate pairwise distance of reads
        if (strcmp(argv[1], "dist")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            UnlabeledValueArg<string> varfileArg("varfile", "path of var file", false, "", "varfile", cmd, true);
            
            cmd.parse(argv2);
            
            cout << "encodefile: " << encodefileArg.getValue() << endl;
            cout << "alignfile: " << alignfileArg.getValue() << endl;
            cout << "outfile: " << outfileArg.getValue() << endl;
            if (varfileArg.getValue()=="")
                cout << "no varfile provided" << endl;
            else
                cout << "varfile: " << varfileArg.getValue() << endl;
            Assembler assembler;
            if (varfileArg.getValue()=="")
                assembler.dist(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue());
            else
                assembler.dist_rdim(encodefileArg.getValue(), alignfileArg.getValue(), varfileArg.getValue(), outfileArg.getValue());
            
        }
        // calculate pairwise jaccard index of reads
        if (strcmp(argv[1], "jaccard")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            ValueArg<double> minjaccardArg("m","minJaccard","minimal jaccard index, default: 0.5", false , 0.5, "minJaccard", cmd);
            
            cmd.parse(argv2);
            
            Assembler assembler;
            assembler.jaccard_index(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(), minjaccardArg.getValue());
        }

        if (strcmp(argv[1], "jaccard_min")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            ValueArg<double> minjaccardArg("m","cutoff","cutoff, default: 0.3", false , 0.3, "cutoff", cmd);
            
            cmd.parse(argv2);
            
            Assembler assembler;
            assembler.jaccard_index_min(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(), minjaccardArg.getValue());
        }

        // pileup encode
        if (strcmp(argv[1], "pileup_var")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            
            int64_t n_reads;
            vector<vector<int> > pu_var = pileup_var(encodefileArg.getValue(), n_reads);
            print_pileup(pu_var, outfileArg.getValue());
        }
        
        // pileup encode
        if (strcmp(argv[1], "pileup_reads")==0) {
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            
            int64_t n_reads;
            vector<vector<int> > pu_reads = pileup_reads(alignfileArg.getValue(), n_reads);
            print_pileup(pu_reads, outfileArg.getValue());
        }
        
        // pileup count
        if (strcmp(argv[1], "pileup_count")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            
            // pileup var
            vector<vector<int> > encode_data;
            loadencodedata(encode_data, encodefileArg.getValue());
            vector<int> pu_var_count = pileup_var_count(encode_data);
            
            // pileup read
            vector<ReadRange> reads_range;
            loadreadsrange(reads_range, alignfileArg.getValue());
            vector<int> pu_reads_count = pileup_reads_m5_count(reads_range);
            
            if (4*(pu_reads_count.size()-1)+3 < pu_var_count.size()-1)
                throw runtime_error("encodefile and alignfile don't match");
            
            // print
            ofstream fs_outfile;
            open_outfile(fs_outfile, outfileArg.getValue());
            for (auto i=0; i < pu_var_count.size(); ++i){
                if (pu_var_count[i]==0)
                    continue;
                double prop = pu_reads_count[int(i/4)] > 0 ? (double)pu_var_count[i] / pu_reads_count[int(i/4)] : 0;
                AlignCoderSNV aligncoder;
                pair<int, char> dc = aligncoder.decode(i);
                fs_outfile << dc.first << '\t' << dc.second << '\t' << i << '\t' << pu_var_count[i] << '\t' << pu_reads_count[int(i/4)] << '\t' << prop << endl;
            }
            fs_outfile.close();
            
        }
        
        // reconstrust reference from m5 file
        if (strcmp(argv[1], "recons_ref") == 0){
            // parse arguments
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            Assembler assembler;
            AlignReaderM5 AlignReaderM5_obj;
            stxxl::vector<Align> align_data;
            AlignReaderM5_obj.read(alignfileArg.getValue(), align_data);
            
            string ref_name; string ref_seq;
            assembler.ref_reconstruct(align_data, ref_name, ref_seq);
            
            ofstream fs_outfile;
            open_outfile(fs_outfile, outfileArg.getValue());
            fs_outfile << ">" << ref_name << endl;
            fs_outfile << ref_seq;
            fs_outfile.close();
        }
        
        // construct haplotype sequence
        if (strcmp(argv[1], "cons_haplo_seq") == 0){
            // parse arguments
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file which encode variants", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file (only 1 chromosome is allowed)", true, "", "reffile", cmd);
            UnlabeledValueArg<string> outdirArg("outdir", "path of output directory", true, "", "outdir", cmd);
            
            cmd.parse(argv2);
            
            // read encodefile
            vector<vector<int> > encode_data;
            loadencodedata(encode_data, encodefileArg.getValue());
            
            // read reference file
            unordered_map<string, string> ref_seq_all;
            read_fasta(reffileArg.getValue(), ref_seq_all);
            if (ref_seq_all.size() != 1)
                throw runtime_error("number of chrosome is not 1.");
            string ref_seq = ref_seq_all.begin()->second; 
            
            // construct haplotype
            system( ("mkdir -p " + outdirArg.getValue()).c_str() );
            
            for (int i = 0; i < (int)encode_data.size(); ++i){
                Assembler assembler;
                string haplo_seq;
                assembler.haplo_seq_construct(encode_data[i], ref_seq, haplo_seq);
                
                string outfile = outdirArg.getValue() + "/haplotype_" + to_string(i) + ".fa";
                ofstream fs_outfile;
                open_outfile(fs_outfile, outfile);
                fs_outfile << ">haplotype_" + to_string(i) << endl;
                fs_outfile << haplo_seq << endl;
                fs_outfile.close();
            }
        }

        if (strcmp(argv[1], "ann")==0){
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> varfileArg("varfile", "path of variant file", true, "", "varfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<int> mincvgArg("c","mincvg","minimal coverage, default: 20", false , 20, "mincvg", cmd);
            ValueArg<double> minpropArg("p","minprop","minimal frequency, default: 0.2", false , 0.2, "minprop", cmd);
            ValueArg<double> maxpropArg("q","maxprop","maximal frequency, default: 0.7", false , 0.7, "maxprop", cmd);
            ValueArg<int> topnArg("t","topn","number of initial neighbors, default: 30", false , 30, "topn", cmd);
            ValueArg<int> maxnnArg("m","maxnn","maximal number of neighbors, default: 200", false , 200, "maxnn", cmd);
            ValueArg<double> maxdistArg("d","maxdist","maximal hamming distance of the initial neighbors, default: 0.02", false , 0.02, "maxprop", cmd);
            
            cmd.parse(argv2);
            cout << "mincvg = " << mincvgArg.getValue() << endl;
            cout << "minprop = " << minpropArg.getValue() << endl;
            cout << "maxprop = " << maxpropArg.getValue() << endl;
            cout << "topn = " << topnArg.getValue() << endl;
            cout << "maxnn = " << maxnnArg.getValue() << endl;
            cout << "maxdist = " << maxdistArg.getValue() << endl;
            
            Assembler assembler;
            assembler.ann_clust(encodefileArg.getValue(), alignfileArg.getValue(), varfileArg.getValue(), mincvgArg.getValue(),
                                minpropArg.getValue(), maxpropArg.getValue(), topnArg.getValue(), maxnnArg.getValue(), maxdistArg.getValue());
            vector<int64_t> idx;
            assembler.find_nccontigs(idx);
            assembler.print_rl_ann_clust(outfileArg.getValue() + ".igda_tmp", false, idx);
            
            string cmd = "sort -u -s -k2n -k3n " + outfileArg.getValue() + ".igda_tmp" + " > " + outfileArg.getValue();
            cout << cmd << endl; system(cmd.c_str());
            
            cmd = "rm -f " + outfileArg.getValue() + ".igda_tmp";
            cout << cmd << endl; system(cmd.c_str());
            
            assembler.print_nc_reads_id(outfileArg.getValue() + ".nc_idx");
            
        }
        
        if (strcmp(argv[1], "samtom5")==0){
            UnlabeledValueArg<string> samfileArg("samfile", "path of sam file", true, "", "samfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> m5fileArg("m5file", "path of m5 file", true, "", "m5file", cmd);
            
            cmd.parse(argv2);
            
            AlignReaderSam alignreadersam;
            alignreadersam.samtom5(samfileArg.getValue(), reffileArg.getValue(), m5fileArg.getValue());
        }
        
        // detect single loci
        if (strcmp(argv[1], "detectsingle")==0){
            UnlabeledValueArg<string> pileupfileArg("pileupfile", "path of pileup file", true, "", "pileupfile", cmd);
            UnlabeledValueArg<string> contextfileArg("contextfile", "path of context file", true, "", "contextfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<double> minbfArg("b","minbf","minimal Bayes factor, default: 50", false , 50, "minbf", cmd);
            ValueArg<double> minfreqArg("f","minfreq","minimal frequency, default: 0.01", false , 0.01, "minfreq", cmd);
            ValueArg<int> mincvgArg("c","mincvg","minimal coverage, default: 20", false , 20, "mincvg", cmd);
            ValueArg<int> mincontextcvgArg("v","mincontextcvg","minimal context total coverage, default: 500", false , 500, "mincontextcvg", cmd);
            
            cmd.parse(argv2);
            DetectSingleSNV detectsinglesnv;
            DetectSingle *p_detectsingle = &detectsinglesnv;
            p_detectsingle->loadcontexteffect(contextfileArg.getValue(), mincontextcvgArg.getValue());
            p_detectsingle->detect(pileupfileArg.getValue(), outfileArg.getValue(), log(minbfArg.getValue()), minfreqArg.getValue(), mincvgArg.getValue());
            
        }
        
        // recode
        if (strcmp(argv[1], "recode")==0){
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of alignment file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> varfileArg("varfile", "path of variants file", true, "", "varfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<int> leftlenArg("l","leftlen","length of the upstream context, default: 10", false , 10, "leftlen", cmd);
            ValueArg<int> rightlenArg("r","rightlen","length of the downstream context, default: 10", false , 10, "rightlen", cmd);
            
            SwitchArg isreportrefArg("f", "ref", "encode reference", cmd, false);
            
            cmd.parse(argv2);
            
            AlignCoderSNV aligncodersnv;
            AlignReaderM5 alignreaderm5;
            AlignCoder *p_aligncoder = &aligncodersnv;
            p_aligncoder->setAlignReader(&alignreaderm5);

            p_aligncoder->recode(alignfileArg.getValue(), varfileArg.getValue(), outfileArg.getValue(), leftlenArg.getValue(), rightlenArg.getValue(), isreportrefArg.getValue());
        }

        
    }
    catch(const std::overflow_error& e) {
        cerr << "overflow_error: " << e.what() << endl;
    }
    catch(const std::runtime_error& e) {
        cerr << "runtime_error: " << e.what() << endl;
    }
    catch(const std::exception& e) {
        cerr << "expection: " << e.what() << endl;
    }
    catch(...) {
        
    }
    return 0;
}

#endif
