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
#include "../src/modules/dforest/dforestsnvstl.h"
#include "../src/modules/errormodel/errormodelsnv.h"
#include "../src/modules/assemble/assembler.h"
#include "../src/modules/detectsingle/detectsinglesnv.h"
#include "../src/modules/rsm/rsmsnv.h"
#include "./misc/misc.h"
#include "./misc/permute_reads.h"
#include "./misc/merge_data.h"
#include "./misc/var2vcf.h"

#ifdef _UNITTEST

//using namespace boost::filesystem;

int main(int argc, char* argv[]) {
    //if (!exists("../results"))
    //    create_directory("../results");
    int result = Catch::Session().run( argc, argv );
    return result;
}

#else

using namespace TCLAP;

void print_usage()
{
    cout << "Version 0.9.3" << endl << endl;
    cout << "To detect low-frequency SNVs, use:" << endl;
    cout << "(PacBio data) igda_pipe_detect_pb infile(bam or sam file) reffile contextmodel outdir" << endl;
    cout << "(Nanopore data) igda_pipe_detect_ont infile(bam or sam file) reffile contextmodel outdir" << endl << endl;
    
    cout << "to detect combination of  low-frequency SNVs, use:" << endl;
    cout << "(PacBio data) igda_pipe_phase_pb indir(output of igda_pipe_detect_pb) reffile outdir" << endl;
    cout << "(Nanopore data) igda_pipe_phase_ont indir(output of igda_pipe_detect_ont) reffile outdir" << endl << endl;
    
    cout << "Please note that reffile is the reference fasta file. Current version assumes there is only one contig in samfile and reffile." << endl << endl;
    cout << "contextmodel is the context effect model trained on independent data. They can be download in https://github.com/zhixingfeng/igda_contextmodel" << endl;
    
    cout << "Output format:" << endl;
    cout << "For detecting low-frequency SNVs, realign.var in outdir is the final result." << endl;
    cout << "Column 1 is the locus of detected SNVs." << endl;
    cout << "Column 2 is the alternative base of detected SNVs." << endl;
    cout << "The other columns are reserved for internal use" << endl << endl;
    
    cout << "For detecting combinations of low-frequency SNVs, realign.ann.tested.ft.count.ft.assembled.count.nc.ft in outdir is the final result." << endl;
    cout << "Each row is a contig" << endl;
    cout << "Column 1 is the SNVs of the contigs. It is encoded, for each integer x, floor(x/4) = 0-based locus, and x modulo 4 = base" << endl;
    cout << "Column 2 is start locus (0-based)" << endl;
    cout << "Column 3 is end locus (0-based)" << endl;
    cout << "Column 4 is number of reads aligned to the contig" << endl;
    cout << "Column 5 is coverage of the contig" << endl;
    cout << "The other columns are reserved for internal use" << endl << endl;
    
    cout << "For any questions, inquiry of source code or example data, contact zhixing.feng@mssm.edu" << endl;
    /*cout << "igda [command]" << endl;
    cout << "command = encode : binary coding SNVs of each aligned read." << endl;
    cout << "          cmpreads : pairwise comparison of encoded reads."<< endl;
    cout << "          dforest : detecting SNVs by dforest algorithm (unable to detect SNVs if no reads cover multiple real SNVs)." << endl;
    cout << "          contexteffect : pileup reads and get context effect." << endl;
    cout << "          detectsingle : detecting SNVs locus by locus (lower sensitivity compared to dforest but able to detect SNVs if no reads cover multiple real SNVs)." << endl;
    cout << "          rdim : removing undetected SNVs from encoded reads." << endl;
    cout << "          ann : clustering reads by adaptive nearest neighbor clustering." << endl;*/
    
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
            
            ValueArg<int> nthreadArg("n","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            
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
                    if (nthreadArg.getValue() == 1){
                        cmpreads_topn(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                                      topnArg.getValue(), overlapArg.getValue(), true, !istextArg.getValue(),
                                      isreadIDArg.getValue(), !isnocondprobArg.getValue());
                    }else{
                        cmpreads_topn_multithread(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(), nthreadArg.getValue(),
                                      topnArg.getValue(), overlapArg.getValue(), true, !istextArg.getValue(),
                                      isreadIDArg.getValue(), !isnocondprobArg.getValue());
                    }
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
        
        // random subspace maximization
        if (strcmp(argv[1], "rsm")==0){
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> cmpreadsfileArg("cmpreadsfile", "path of cmpreads file", true, "", "cmpreadsfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<int> minreadsArg("r","minreads","minimal number of reads in a node, default: 12", false , 12, "minreads", cmd);
            ValueArg<int> maxdepthArg("d","maxdepth","maximal depth of a tree, default: 1000", false , 1000, "maxdepth", cmd);
            ValueArg<double> minfreqArg("f","minfreq","minimal frequency: 0.0", false , 0.0, "minfreq", cmd);
            ValueArg<double> maxfreqArg("q","maxfreq","maximal frequency: 1.0 (substantially increase speed if it is small, but restrict the maximal conditional frequency)", false , 1.0, "maxfreq", cmd);
            ValueArg<int> minhomoArg("m","minhomo","minimal homopolymer blocks distance between linked loci, default: 15", false , 15, "minhomo", cmd);
            ValueArg<int> nthreadArg("n","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            
            SwitchArg isinterArg("i", "intermediate", "output intermediate results", cmd, false);
            
            cmd.parse(argv2);
            cout << "minreads = " << minreadsArg.getValue() << endl;
            cout << "maxdepth = " << maxdepthArg.getValue() << endl;
            cout << "minfreq = " << minfreqArg.getValue() << endl;
            cout << "maxfreq = " << maxfreqArg.getValue() << endl;
            cout << "minhomo = " << minhomoArg.getValue() << endl;
            cout << "nthread = " << nthreadArg.getValue() << endl;
            
            if (minfreqArg.getValue() > maxfreqArg.getValue())
                throw runtime_error("minfreqArg.getValue() > maxfreqArg.getValue()");
            
            AlignReaderM5 alignreaderm5;
            AlignCoderSNV aligncoder;
            
            RSMsnv rsmsnv(&alignreaderm5, &aligncoder);
            
            rsmsnv.load_homo_blocks(reffileArg.getValue());
            rsmsnv.run(encodefileArg.getValue(), encodefileArg.getValue()+".ref", alignfileArg.getValue(), cmpreadsfileArg.getValue(),
                       outfileArg.getValue(), minreadsArg.getValue(), maxdepthArg.getValue(), nthreadArg.getValue(), minfreqArg.getValue(),
                       maxfreqArg.getValue(), minhomoArg.getValue(), isinterArg.getValue());
        }
        
        // dforest algorithm
        if (strcmp(argv[1], "dforest")==0){
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> cmpreadsfileArg("cmpreadsfile", "path of cmpreads file", true, "", "cmpreadsfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            UnlabeledValueArg<string> tmpdirArg("tmpdir", "temporary directory", true, "", "tmpdir", cmd);
            
            ValueArg<int> minreadsArg("r","minreads","minimal number of reads in a node, default: 25", false , 25, "minreads", cmd);
            ValueArg<int> maxdepthArg("d","maxdepth","maximal depth of a tree, default: 1000", false , 1000, "maxdepth", cmd);
            ValueArg<double> minfreqArg("f","minfreq","minimal frequency: 0.0", false , 0.0, "minfreq", cmd);
            ValueArg<double> maxfreqArg("q","maxfreq","maximal frequency: 1.0 (substantially increase speed if it is small, but restrict the maximal conditional frequency)", false , 1.0, "maxfreq", cmd);
            ValueArg<int> minhomoArg("m","minhomo","minimal homopolymer blocks distance between linked loci, default: 15", false , 15, "minhomo", cmd);
            ValueArg<int> nthreadArg("n","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            ValueArg<string> varfileArg("v","varfile","variant file", false , "", "varfile", cmd);
            
            SwitchArg islegacyArg("l", "legacy", "use the legacy algorithm (no stxxl) to run dforest", cmd, false);
            SwitchArg isinterArg("i", "intermediate", "output intermediate results", cmd, false);
            
            cmd.parse(argv2);
            cout << "minreads = " << minreadsArg.getValue() << endl;
            cout << "maxdepth = " << maxdepthArg.getValue() << endl;
            cout << "minfreq = " << minfreqArg.getValue() << endl;
            cout << "maxfreq = " << maxfreqArg.getValue() << endl;
            cout << "minhomo = " << minhomoArg.getValue() << endl;
            cout << "nthread = " << nthreadArg.getValue() << endl;
            
            if (minfreqArg.getValue() > maxfreqArg.getValue())
                throw runtime_error("minfreqArg.getValue() > maxfreqArg.getValue()");
            
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
            //DForestSNVMax forestsnvmax(p_alignreader, &aligncoder);
            DForestSNVSTL forestsnvstxxl(p_alignreader, &aligncoder);
            
            DForest *ptr_forest;
            ptr_forest = &forestsnvstxxl;
            
            string shell_cmd = "mkdir -p " + tmpdirArg.getValue();
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            
            ptr_forest->load_homo_blocks(reffileArg.getValue());
            ptr_forest->run(encodefileArg.getValue(), alignfileArg.getValue(), cmpreadsfileArg.getValue(), outfileArg.getValue(), tmpdirArg.getValue(),
                            minreadsArg.getValue(), maxdepthArg.getValue(), nthreadArg.getValue(), minfreqArg.getValue(), maxfreqArg.getValue(), minhomoArg.getValue(),
                            isinterArg.getValue(), varfileArg.getValue());
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
            UnlabeledValueArg<string> infileArg("infile", "path of input file (m5 file if no -p, pileup file if -p)", true, "", "infile", cmd);
            UnlabeledValueArg<string> outprefixArg("outprefix", "prefix of output files", true, "", "outprefix", cmd);
            
            ValueArg<int> leftlenArg("l","leftlen","length of the upstream context, default: 1", false , 1, "leftlen", cmd);
            ValueArg<int> rightlenArg("r","rightlen","length of the downstream context, default: 1", false , 1, "rightlen", cmd);
            
            SwitchArg ispileupArg("p", "ispileup", "is use pileup file as input", cmd, false);
            
            cmd.parse(argv2);
            ErrorModelSNV errormodel;
            errormodel.set_context_size(leftlenArg.getValue(), rightlenArg.getValue());
            
            if (ispileupArg.getValue()){
                vector<BaseFreq> pileup = errormodel.load_pileup(infileArg.getValue());
                ContextEffect context_effect;
                errormodel.get_context_effect(pileup, context_effect);
                errormodel.print_context_effect(outprefixArg.getValue() + ".context", context_effect);
            }else{
                errormodel.learn(infileArg.getValue(), outprefixArg.getValue());
            }
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

        // pileup qv
        if (strcmp(argv[1], "pileup_qv")==0) {
            UnlabeledValueArg<string> samfileArg("samfile", "path of SAM file", true, "", "samfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference fasta file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            SwitchArg iscountArg("c", "iscount", "is use average qv instead of pileup all qv", cmd, false);
            SwitchArg isvarArg("v", "isvar", "is variants only", cmd, false);
            
            cmd.parse(argv2);
            
            if (iscountArg.getValue()){
                vector<vector<double> > pu_qv_count = pileup_qv_count(samfileArg.getValue(), reffileArg.getValue(), isvarArg.getValue());
                print_pileup_qv_count(pu_qv_count, outfileArg.getValue());
            }else{
                vector<vector<pair<int64_t, double> > > pu_qv = pileup_qv(samfileArg.getValue(), reffileArg.getValue(), isvarArg.getValue());
                print_pileup_qv(pu_qv, outfileArg.getValue());
            }
            
        }
        
        // pileup
        if (strcmp(argv[1], "pileup")==0) {
            UnlabeledValueArg<string> m5fileArg("m5file", "path of m5 file", true, "", "m5file", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            
            ErrorModelSNV errormodel;
            int g_size = errormodel.get_genomesize(m5fileArg.getValue());
            vector<BaseFreq> pileup(g_size, BaseFreq());
            errormodel.pileup_reads(m5fileArg.getValue(), pileup);
            errormodel.print_pileup(outfileArg.getValue(), pileup);
           
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
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
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
        
        // pileup_count_to_context
        if (strcmp(argv[1], "pileup_count_to_context")==0) {
            UnlabeledValueArg<string> pu_count_fileArg("pu_count_file", "path of pu_count_file", true, "", "pu_count_file", cmd);
            UnlabeledValueArg<string> pu_fileArg("pu_file", "path of pu_file (the .pileup file in output of igda contexteffect)", true, "", "pu_file", cmd);
            UnlabeledValueArg<string> out_fileArg("out_file", "path of out_file", true, "", "out_file", cmd);
            
            cmd.parse(argv2);
            
            ErrorModelSNV model;
            model.pileup_count_to_context(pu_count_fileArg.getValue(), pu_fileArg.getValue(), out_fileArg.getValue());
                        
        }

        // reconstrust reference from m5 file
        if (strcmp(argv[1], "recons_ref") == 0){
            // parse arguments
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            Assembler assembler;
            AlignReaderM5 AlignReaderM5_obj;
            vector<Align> align_data;
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

        if (strcmp(argv[1], "find_nccontigs")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            ValueArg<double> minpropArg("p","minprop","minimal proportion of overlaped SNVs (default 0.5)", false , 0.5, "minprop", cmd);
            ValueArg<double> minjaccardpArg("j","minjaccard","minimal jaccard index of the two contigs, default: 2 (i.e. not use jaccard)", false , 2, "minjaccard", cmd);
            
            cmd.parse(argv2);
            Assembler assembler;
            assembler.read_ann_results(annfileArg.getValue());
            vector<int64_t> idx;
            assembler.find_nccontigs(idx, minpropArg.getValue(), minjaccardpArg.getValue());
            assembler.print_rl_ann_clust(outfileArg.getValue() + ".igda_tmp", true, idx);
            string cmd = "sort -u -s -k2,2n -k3,3n -k1,1 " + outfileArg.getValue() + ".igda_tmp" + " > " + outfileArg.getValue();
            cout << cmd << endl; system(cmd.c_str());
            cmd = "rm -f " + outfileArg.getValue() + ".igda_tmp";
            cout << cmd << endl; system(cmd.c_str());
        }
        if (strcmp(argv[1], "ann")==0){
            UnlabeledValueArg<string> recodefileArg("recodefile", "path of recode file", true, "", "recodefile", cmd);
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> varfileArg("varfile", "path of variant file", true, "", "varfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<int> mincvgArg("c","mincvg","minimal coverage, default: 10", false , 10, "mincvg", cmd);
            ValueArg<double> minpropArg("p","minprop","minimal frequency, default: 0.2", false , 0.2, "minprop", cmd);
            ValueArg<double> maxpropArg("q","maxprop","maximal frequency, default: 0.8", false , 0.8, "maxprop", cmd);
            ValueArg<int> topnArg("t","topn","number of initial neighbors, default: 12", false , 25, "topn", cmd);
            ValueArg<int> maxnnArg("m","maxnn","maximal number of neighbors, default: 30", false , 50, "maxnn", cmd);
            ValueArg<double> minjaccardArg("j","minjaccard","minimal jaccard index of the initial neighbors, default: 0.5", false , 0.5, "minjaccard", cmd);
            ValueArg<int> maxiterArg("b","b","maximal number of iterations, default: 1", false , 1, "b", cmd);
            ValueArg<int> nthreadArg("n","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            
            SwitchArg iscorrectArg("r", "correct", "is use corrected reads as seeds", cmd, false);
            SwitchArg islegacyArg("l", "legacy", "is use legacy version (no recoding)", cmd, false);
            SwitchArg ishangArg("g", "ishang", "is only use reads whose start is smaller than the seed as neighbors", cmd, false);
            SwitchArg isrecodeArg("e", "isrecode", "is use recode_data instead of encode_data to calculate jaccard index", cmd, false);
            
            cmd.parse(argv2);
            cout << "mincvg = " << mincvgArg.getValue() << endl;
            cout << "minprop = " << minpropArg.getValue() << endl;
            cout << "maxprop = " << maxpropArg.getValue() << endl;
            cout << "topn = " << topnArg.getValue() << endl;
            cout << "maxnn = " << maxnnArg.getValue() << endl;
            cout << "minjaccard = " << minjaccardArg.getValue() << endl;
            
            Assembler assembler;

            if (islegacyArg.getValue()){
                assembler.ann_clust(encodefileArg.getValue(), alignfileArg.getValue(), varfileArg.getValue(), mincvgArg.getValue(),
                                minpropArg.getValue(), maxpropArg.getValue(), topnArg.getValue(), maxnnArg.getValue(), minjaccardArg.getValue());
            }else{
                if (nthreadArg.getValue() == 1){
                    assembler.ann_clust_recode(recodefileArg.getValue(), recodefileArg.getValue() + ".ref", encodefileArg.getValue(), alignfileArg.getValue(), varfileArg.getValue(), mincvgArg.getValue(),
                                               minpropArg.getValue(), maxpropArg.getValue(), topnArg.getValue(), maxnnArg.getValue(), minjaccardArg.getValue(), iscorrectArg.getValue(), true, maxiterArg.getValue(), isrecodeArg.getValue());
                }else{
                    assembler.ann_clust_recode_multithread(recodefileArg.getValue(), recodefileArg.getValue() + ".ref", encodefileArg.getValue(), alignfileArg.getValue(), mincvgArg.getValue(),
                                               minpropArg.getValue(), maxpropArg.getValue(), topnArg.getValue(), maxnnArg.getValue(), minjaccardArg.getValue(), iscorrectArg.getValue(), true, maxiterArg.getValue(), isrecodeArg.getValue(), nthreadArg.getValue());
                }
            }
            vector<int64_t> idx;
            assembler.print_rl_ann_clust(outfileArg.getValue() + ".raw", true);
            assembler.find_nccontigs(idx);
            assembler.print_rl_ann_clust(outfileArg.getValue() + ".igda_tmp", true, idx);
            string cmd = "sort -u -s -k2,2n -k3,3n -k1,1 " + outfileArg.getValue() + ".igda_tmp" + " > " + outfileArg.getValue();
            cout << cmd << endl; system(cmd.c_str());
            
            cmd = "rm -f " + outfileArg.getValue() + ".igda_tmp";
            cout << cmd << endl; system(cmd.c_str());
            
        }
        
        if (strcmp(argv[1], "test_contigs")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            
            ValueArg<double> minlogbfArg("b","minlogbf","minimal logrithm of bayes factors, default: 50", false , 50, "minlogbf", cmd);
            ValueArg<double> maxlociArg("l","maxloci","maximal number of loci to avoid testing, default: 10", false , 10, "maxloci", cmd);
            ValueArg<double> minrrArg("r","minrr","minimal relative risk, default: 10", false , 10, "minrr", cmd);
            
            ValueArg<double> alphaArg("a","alpha","alpha of null beta distribution, default: 1.332824", false , 1.332824, "alpha", cmd);
            ValueArg<double> betaArg("t","beta","beta of null beta distribution, default: 89.04769", false , 89.04769, "beta", cmd);
            
            
            cmd.parse(argv2);
            Assembler assembler;
            
            cout << "load ann" << endl;
            assembler.read_ann_results(annfileArg.getValue());

            cout << "load recode_data" << endl;
            vector<vector<int> > recode_data;
            loadencodedata(recode_data, encodefileArg.getValue());
            
            cout << "load recode_ref_data" << endl;
            vector<vector<int> > recode_ref_data;
            loadencodedata(recode_ref_data, encodefileArg.getValue() + ".ref");
            
            cout << "load reads_range" << endl;
            vector<ReadRange> reads_range;
            loadreadsrange(reads_range, alignfileArg.getValue());
            
            cout << "load ref_file" << endl;
            assembler.load_homo_blocks(reffileArg.getValue());
            
            cout << "test contigs" << endl;
            assembler.set_null_betadist(alphaArg.getValue(), betaArg.getValue());
            assembler.test_contigs(recode_data, recode_ref_data, reads_range);
            assembler.print_rl_ann_clust(annfileArg.getValue() + ".tested", true);
            
            cout << "filter contigs" << endl;
            assembler.filter_ann(annfileArg.getValue() + ".tested", minlogbfArg.getValue(), maxlociArg.getValue(), minrrArg.getValue());
            
        }
        
        if (strcmp(argv[1], "test_contigs_pairwise")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> recodefileArg("recodefile", "path of recode file", true, "", "recodefile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
           
            ValueArg<double> alphaArg("a","alpha","alpha of null beta distribution, default: 1.332824", false , 1.332824, "alpha", cmd);
            ValueArg<double> betaArg("t","beta","beta of null beta distribution, default: 89.04769", false , 89.04769, "beta", cmd);
            
            ValueArg<double> minlogbfArg("b","minlogbf","minimal logrithm of bayes factors, default: 50", false , 50, "minlogbf", cmd);
            ValueArg<double> maxlociArg("l","maxloci","maximal number of loci to avoid testing, default: 3", false , 3, "maxloci", cmd);
            ValueArg<double> mincvgArg("c","mincvg","minimal coverage, default: 10", false , 10, "mincvg", cmd);
            ValueArg<double> minrrArg("r","minrr","minimal relative risk, default: 10", false , 10, "minrr", cmd);
            
            
            cmd.parse(argv2);
         
            Assembler assembler;
            assembler.load_homo_blocks(reffileArg.getValue());
            assembler.set_null_betadist(alphaArg.getValue(), betaArg.getValue());
            assembler.test_contigs_pairwise(annfileArg.getValue(), recodefileArg.getValue(), annfileArg.getValue() + ".ft",
                                            minlogbfArg.getValue(), maxlociArg.getValue(), mincvgArg.getValue(), minrrArg.getValue());
            
        }
        
        if (strcmp(argv[1], "polish")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            UnlabeledValueArg<string> tmpdirArg("tmpdir", "temporary directory", true, "", "tmpdir", cmd);
            
            ValueArg<int> minreadsArg("r","minreads","minimal number of reads in a node, default: 10", false , 10, "minreads", cmd);
            ValueArg<double> minfreqArg("f","minfreq","minimal frequency: 0.3", false , 0.3, "minfreq", cmd);
            ValueArg<int> minhomoArg("m","minhomo","minimal homopolymer blocks distance between linked loci, default: 1", false , 1, "minhomo", cmd);
            ValueArg<string> varfileArg("v","varfile","variant file", false , "", "varfile", cmd);
            
            cmd.parse(argv2);
            cout << "minreads = " << minreadsArg.getValue() << endl;
            cout << "minfreq = " << minfreqArg.getValue() << endl;
            cout << "minhomo = " << minhomoArg.getValue() << endl;
            
            Assembler assembler;
            assembler.polish(annfileArg.getValue(), encodefileArg.getValue(), alignfileArg.getValue(),
                             reffileArg.getValue(), outfileArg.getValue(), tmpdirArg.getValue(),
                             minfreqArg.getValue(), minreadsArg.getValue(), minhomoArg.getValue(), varfileArg.getValue());
           
        }
        
        if (strcmp(argv[1], "correct_contigs")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<double> minoverlapArg("l","minoverlap","minimal overlap ratio between two contigs, default: 0", false , 0, "minoverlap", cmd);
            
            
            cmd.parse(argv2);
            Assembler assembler;
            assembler.correct_contigs(annfileArg.getValue(), outfileArg.getValue(), minoverlapArg.getValue());
           
            
        }
        
        
        if (strcmp(argv[1], "abundance")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> recodefileArg("recodefile", "path of recode file", true, "", "recodefile", cmd);
            UnlabeledValueArg<string> m5fileArg("m5file", "path of m5 file", true, "", "m5file", cmd);
            
            cmd.parse(argv2);
            
            Assembler assembler;
            cout << "load ann" << endl;
            assembler.read_ann_results(annfileArg.getValue());
            
            cout << "load recode_data" << endl;
            vector<vector<int> > recode_data;
            loadencodedata(recode_data, recodefileArg.getValue());
            
            cout << "load reads_range" << endl;
            vector<ReadRange> reads_range;
            loadreadsrange(reads_range, m5fileArg.getValue());
            
            assembler.assign_reads_to_contigs(recode_data, reads_range);
            
            assembler.print_rl_ann_clust(annfileArg.getValue() + ".count", true);
          
        }
        
        if (strcmp(argv[1], "tred")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            
            ValueArg<double> minpropArg("p","minprop","minimal proportion of common variants between two contigs, default: 0.5", false , 0.5, "minprop", cmd);
            ValueArg<double> minlenpropArg("l","minlenprop","minimal proportion of overlaping length between two contigs, default: 0.5", false , 0.5, "minlenprop", cmd);
            ValueArg<double> minjaccardpArg("j","minjaccard","minimal jaccard index of the two contigs, default: 2 (i.e. not use jaccard)", false , 2, "minjaccard", cmd);
            
            SwitchArg isorfArg("r", "or", "is use or in overlap cutoff", cmd, false);
            
            cmd.parse(argv2);
            
            // get overlap graph
            Assembler assembler;
            Graph gp_bgl;
            assembler.ann_to_graph(gp_bgl, annfileArg.getValue(), minpropArg.getValue(), minlenpropArg.getValue(), minjaccardpArg.getValue(), isorfArg.getValue());
            
            ofstream fs_graph(annfileArg.getValue() + ".dot");
            boost::write_graphviz(fs_graph, gp_bgl);
            fs_graph.close();
            
            // transitive reduction
            IGDA_Graph gp;
            load_igda_graph_from_file(gp, annfileArg.getValue() + ".dot", annfileArg.getValue());
            
            IGDA_Graph gp_tred;
            igda_tred(gp, gp_tred);
            
            save_igda_graph_to_file(gp_tred, annfileArg.getValue() + ".tred.dot");
            
            /*Assembler assembler;
            Graph gp;
            assembler.ann_to_graph(gp, annfileArg.getValue(), minpropArg.getValue(), minlenpropArg.getValue());
            
            ofstream fs_graph(annfileArg.getValue() + ".dot");
            boost::write_graphviz(fs_graph, gp);
            fs_graph.cloase();
            
            string cmd = "tred " + annfileArg.getValue() + ".dot";
            cmd = cmd + " > " + annfileArg.getValue() + ".tred.dot";
            cout << cmd << endl;
            system(cmd.c_str());*/
        }
        
        if (strcmp(argv[1], "tred_legacy")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            
            ValueArg<double> minpropArg("p","minprop","minimal proportion of common variants between two contigs, default: 0.5", false , 0.5, "minprop", cmd);
            ValueArg<double> minlenpropArg("l","minlenprop","minimal proportion of overlaping length between two contigs, default: 0.5", false , 0.5, "minlenprop", cmd);
            
            cmd.parse(argv2);
            
            Assembler assembler;
            Graph gp;
            assembler.ann_to_graph(gp, annfileArg.getValue(), minpropArg.getValue(), minlenpropArg.getValue());
            
            ofstream fs_graph(annfileArg.getValue() + ".dot");
            boost::write_graphviz(fs_graph, gp);
            fs_graph.close();
            
            string cmd = "tred " + annfileArg.getValue() + ".dot";
            cmd = cmd + " > " + annfileArg.getValue() + ".tred.dot";
            cout << cmd << endl;
            system(cmd.c_str());
        }
        
        if (strcmp(argv[1], "assemble")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> dotfileArg("dotfile", "path of dot file (with transitive reduction)", true, "", "dotfile", cmd);
            
            SwitchArg legacyfArg("l", "legacy", "is use legacy algorithm", cmd, false);
            
            cmd.parse(argv2);
            
            if (legacyfArg.getValue()){
                cout << "use the legacy assembly algorithm" << endl;
                Assembler assembler;
                assembler.read_ann_results(annfileArg.getValue());
                Graph gp;
                read_dot_file(gp, dotfileArg.getValue());
                assembler.assemble(gp, annfileArg.getValue() + ".assembled");
            }else{
                Assembler assembler;
                assembler.read_ann_results(annfileArg.getValue());
                IGDA_Graph gp;
                load_igda_graph_from_file(gp, dotfileArg.getValue(), annfileArg.getValue());
                assembler.assemble_unambiguous(gp, annfileArg.getValue() + ".assembled");
            }
            
        }

        if (strcmp(argv[1], "samtom5")==0){
            UnlabeledValueArg<string> samfileArg("samfile", "path of sam file", true, "", "samfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> m5fileArg("m5file", "path of m5 file", true, "", "m5file", cmd);
            
            ValueArg<int> minlenArg("l","minlen","minimal length of alignment", false , 1000, "minlen", cmd);
            
            cmd.parse(argv2);
            
            AlignReaderSam alignreadersam;
            alignreadersam.samtom5(samfileArg.getValue(), reffileArg.getValue(), m5fileArg.getValue(), minlenArg.getValue());
        }
        
        if (strcmp(argv[1], "samtom5qv")==0){
            UnlabeledValueArg<string> samfileArg("samfile", "path of sam file", true, "", "samfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> m5qvfileArg("m5qvfile", "path of m5qv file", true, "", "m5qvfile", cmd);
            
            ValueArg<int> minlenArg("l","minlen","minimal length of alignment", false , 1000, "minlen", cmd);
            
            cmd.parse(argv2);
            
            AlignReaderSam alignreadersam;
            alignreadersam.samtom5qv(samfileArg.getValue(), reffileArg.getValue(), m5qvfileArg.getValue(), minlenArg.getValue());
        }
        
        if (strcmp(argv[1], "getbamchrrange")==0){
            UnlabeledValueArg<string> samfileArg("samfile", "path of sam(bam) file", true, "", "samfile", cmd);
            UnlabeledValueArg<string> reffileArg("reffile", "path of reference file", true, "", "reffile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            
            AlignReaderSam alignreadersam;
            alignreadersam.getchrrange(samfileArg.getValue(), reffileArg.getValue(), outfileArg.getValue());
        }
        
        // detect single loci
        if (strcmp(argv[1], "detectsingle")==0){
            UnlabeledValueArg<string> pileupfileArg("pileupfile", "path of pileup file", true, "", "pileupfile", cmd);
            UnlabeledValueArg<string> contextfileArg("contextfile", "path of context file", true, "", "contextfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<double> minbfArg("b","minbf","minimal log Bayes factor, default: 10", false , 10, "minbf", cmd);
            ValueArg<double> minfreqArg("f","minfreq","minimal frequency, default: 0.01", false , 0.01, "minfreq", cmd);
            ValueArg<int> mincvgArg("c","mincvg","minimal coverage, default: 20", false , 15, "mincvg", cmd);
            ValueArg<int> mincontextcvgArg("v","mincontextcvg","minimal context total coverage, default: 0", false , 0, "mincontextcvg", cmd);
            
            cmd.parse(argv2);
            DetectSingleSNV detectsinglesnv;
            DetectSingle *p_detectsingle = &detectsinglesnv;
            p_detectsingle->loadcontexteffect(contextfileArg.getValue(), mincontextcvgArg.getValue());
            p_detectsingle->detect(pileupfileArg.getValue(), outfileArg.getValue(), minbfArg.getValue(), minfreqArg.getValue(), mincvgArg.getValue());
            
        }
        
        // recode
        if (strcmp(argv[1], "recode")==0){
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of alignment file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> varfileArg("varfile", "path of variants file", true, "", "varfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<int> leftlenArg("l","leftlen","length of the upstream context, default: 10", false , 10, "leftlen", cmd);
            ValueArg<int> rightlenArg("r","rightlen","length of the downstream context, default: 10", false , 10, "rightlen", cmd);
            ValueArg<int> nthreadArg("t","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            
            SwitchArg isnorefArg("n", "noref", "is not encode reference", cmd, false);
            //SwitchArg islegacyArg("y", "legacy", "is use legacy algorithm", cmd, false);
            
            cmd.parse(argv2);
            
            AlignCoderSNV aligncodersnv;
            AlignReaderM5 alignreaderm5;
            AlignCoder *p_aligncoder = &aligncodersnv;
            
            if (nthreadArg.getValue() <= 1){
                p_aligncoder->recode(alignfileArg.getValue(), varfileArg.getValue(), outfileArg.getValue(), leftlenArg.getValue(), rightlenArg.getValue(), !isnorefArg.getValue());
            }else{
                p_aligncoder->recode_multithread(alignfileArg.getValue(), varfileArg.getValue(), outfileArg.getValue(), leftlenArg.getValue(), rightlenArg.getValue(), nthreadArg.getValue(),  !isnorefArg.getValue());
            }
            
        }
        
        // permute reads
        if (strcmp(argv[1], "permute")==0){
            UnlabeledValueArg<string> m5fileArg("m5file", "path of m5 file", true, "", "m5file", cmd);
            UnlabeledValueArg<string> pufileArg("pufile", "path of pileup file", true, "", "pufile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<int> seedArg("s","seed","seed of random number generator, default: 18473", false , 18473, "seed", cmd);
            
            cmd.parse(argv2);
            permute_encodefile(m5fileArg.getValue(), pufileArg.getValue(), outfileArg.getValue(), seedArg.getValue());
            
         }
        
        // get_condprob_threshold
        if (strcmp(argv[1], "get_condprob_threshold")==0){
            UnlabeledValueArg<string> permutefileArg("permuted_file", "path of permuted file", true, "", "permuted_file", cmd);
            UnlabeledValueArg<string> pufileArg("pufile", "path of pileup file", true, "", "pufile", cmd);
           
            ValueArg<double> minprobArg("p","minprob","minimal probablity to consider", false , 0.2, "minprob", cmd);
           
            cmd.parse(argv2);
            get_condprob_threshold(permutefileArg.getValue(), pufileArg.getValue(), minprobArg.getValue());
            
        }
        
        // merge encode
        if (strcmp(argv[1], "merge_encode")==0){
            UnlabeledValueArg<string> m5fofnfileArg("m5_fofn_file", "path of m5 fofn file", true, "", "m5_fofn_file", cmd);
            UnlabeledValueArg<string> encodefofnfileArg("encode_fofn_file", "path of encode fofn file", true, "", "encode_fofn_file", cmd);
            UnlabeledValueArg<string> outfileArg("out_file", "path of output encode file", true, "", "out_file", cmd);
           
            cmd.parse(argv2);
            merge_encode(m5fofnfileArg.getValue(), encodefofnfileArg.getValue(), outfileArg.getValue());
        }
        
        // merge m5
        if (strcmp(argv[1], "merge_m5")==0){
            UnlabeledValueArg<string> m5fofnfileArg("m5_fofn_file", "path of m5 fofn file", true, "", "m5_fofn_file", cmd);
            UnlabeledValueArg<string> readnamefileArg("readname_file", "path of readname list file", true, "", "readname_file", cmd);
            UnlabeledValueArg<string> outfileArg("out_file", "path of output m5 file", true, "", "out_file", cmd);
           
            cmd.parse(argv2);
            merge_m5(m5fofnfileArg.getValue(), readnamefileArg.getValue(), outfileArg.getValue());
        }
        
        // split file
        if (strcmp(argv[1], "split_file")==0){
            UnlabeledValueArg<string> infileArg("infile", "path of input file", true, "", "infile", cmd);
            UnlabeledValueArg<string> outprefixArg("outprefix", "path of output prefix", true, "", "outprefix", cmd);
            
            ValueArg<int> nfileArg("n","nfile","number of files, default: 2", false , 2, "nfile", cmd);
            
            cmd.parse(argv2);
            size_t nlines = get_file_nlines(infileArg.getValue());
            size_t nlines_per_file = ceil ( (double) nlines / nfileArg.getValue());
            split_file(infileArg.getValue(), outprefixArg.getValue(), nlines_per_file);
        }
        
        if (strcmp(argv[1], "accessible_vertices")==0){
            UnlabeledValueArg<string> annfileArg("annfile", "path of ann file", true, "", "annfile", cmd);
            UnlabeledValueArg<string> dotfileArg("dotfile", "path of dot file (with transitive reduction)", true, "", "dotfile", cmd);
            
            ValueArg<int> vertexArg("v","vertex","all connected vertices from the starting vertex", true , 0, "vertex", cmd);
            
            cmd.parse(argv2);
            
            IGDA_Graph gp;
            load_igda_graph_from_file(gp, dotfileArg.getValue(), annfileArg.getValue());
            unordered_set<int64_t> accessible_vertices;
            get_accessible_vertices(gp, accessible_vertices, vertexArg.getValue());
            for (auto v : accessible_vertices)
                cout << v << ' ';
            cout << endl;
            
        }
        
        if (strcmp(argv[1], "getdepth_from_paf")==0){
            UnlabeledValueArg<string> paf_fileArg("paf_file", "path of PAF file", true, "", "paf_file", cmd);
            UnlabeledValueArg<string> cvg_fileArg("cvg_file", "path of output coverage file)", true, "", "cvg_file", cmd);
            
            cmd.parse(argv2);
            
            getdepth_from_paf(paf_fileArg.getValue(), cvg_fileArg.getValue());
        }
        
        if (strcmp(argv[1], "var2vcf")==0){
            UnlabeledValueArg<string> var_fileArg("var_file", "path of var file", true, "", "var_file", cmd);
            UnlabeledValueArg<string> ref_fileArg("ref_file", "path of ref file", true, "", "ref_file", cmd);
            UnlabeledValueArg<string> vcf_fileArg("vcf_file", "path of vcf file", true, "", "vcf_file", cmd);
            
            cmd.parse(argv2);
            var2vcf(var_fileArg.getValue(), ref_fileArg.getValue(), vcf_fileArg.getValue());
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
