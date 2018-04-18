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
#include "../src/modules/dforest/dforestsnv.h"
#include "../src/modules/dforest/dforestsnvfast.h"
#include "../src/modules/dforest/dforestsnvmax.h"
#include "../src/modules/errormodel/errormodelsnv.h"
#include "../src/modules/hclust/hclust.h"
#include "../src/modules/sclust/sclust.h"
#include "../src/modules/assemble/assembler.h"
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
    cout << "command = samtofa, bamtofa, m5tofa, encode, cmpreads, sclust, eval, bin2txt, txt2bin, dforest, sort, filter, contexteffect, merge, mergeall, mask, dist, pileup_var, pileup_reads" << endl;
    cout << "bamtofa: convert bam file to fasta file, convert sequence mapped to negative strand to its reverse complementary sequence" << endl;
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
            ValueArg<int> mArg("m","method","method for encoding. 0: full, 1: SNV", false , 0, "method", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of alignment file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            cmd.parse(argv2);
            
            // set alignreader
            AlignReaderM5 alignreaderm5;
            
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
            
            p_aligncoder->setAlignReader(&alignreaderm5);
            p_aligncoder->encode(alignfileArg.getValue(), outfileArg.getValue());
            
            return 0;
        }
        
        // pairwise compare reads
        if (strcmp(argv[1], "cmpreads")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            ValueArg<int> topnArg("p","topn","select top n candidates of each reads, default: 10", false , 0, "topn", cmd);
            ValueArg<double> overlapArg("l","overlap","minimal overlap of reads, default: 0.25", false , 0.25, "overlap", cmd);
            SwitchArg istextArg("t", "text", "is output text file", cmd, false);
            SwitchArg isreadIDArg("r", "readID", "is print reads ID", cmd, false);
            SwitchArg isnocondprobArg("c", "condprob", "not use conditional probability", cmd, false);
            SwitchArg isdiffArg("d", "diff", "is calculate difference between reads", cmd, false);
            
            cmd.parse(argv2);
            
            if (isdiffArg.getValue()){
                cmpreads_topn_diff(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                                   topnArg.getValue(), overlapArg.getValue());
            }else{
                cmpreads_topn(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                              topnArg.getValue(), overlapArg.getValue(), true, !istextArg.getValue(),
                              isreadIDArg.getValue(), !isnocondprobArg.getValue());
            }
            
            
        }
        
        // split cmpreads_file
        if (strcmp(argv[1], "split") == 0){
            UnlabeledValueArg<string> cmpreadsfileArg("cmpreadsfile", "path of cmpreads file", true, "", "cmpreadsfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            ValueArg<int> maxsubdimArg("s","maxsubdim","maximal subspace dimension, default: 15", false , 15, "maxsubdim", cmd);
            
            cmd.parse(argv2);
            
            SClust sclust;
            sclust.split_subspace(cmpreadsfileArg.getValue(), outfileArg.getValue(), maxsubdimArg.getValue());
        }

        
        // subspace clustering
        if (strcmp(argv[1], "sclust") == 0){
            // parse arguments
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> cmpreadsfileArg("cmpreadsfile", "path of cmpreads file", true, "", "cmpreadsfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            UnlabeledValueArg<string> tmpdirArg("tmpdir", "temporary directory", true, "", "tmpdir", cmd);
            
            ValueArg<int> maxsubdimArg("s","maxsubdim","maximal subspace dimension, default: 5", false , 5, "maxsubdim", cmd);
            ValueArg<double> mincondprobrg("p","minCondProb","minimal conditional probability, default: 0.75", false , 0.75, "minCondProb", cmd);
            ValueArg<int> mincountArg("c","mincount","minimal count of variants: 10", false , 10, "mincount", cmd);
            ValueArg<int> mincvgArg("v","cvg","minimal coverage of subspace: 10", false , 10, "mincvg", cmd);
            ValueArg<int> nthreadArg("n","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            
            cmd.parse(argv2);
            
            // make temporary directory 
            string shell_cmd = "mkdir -p " + tmpdirArg.getValue();
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            
            
            SClust sclust;
            sclust.run(encodefileArg.getValue(), alignfileArg.getValue(), cmpreadsfileArg.getValue(),
                       outfileArg.getValue(), tmpdirArg.getValue(), maxsubdimArg.getValue(),
                       mincondprobrg.getValue(), mincountArg.getValue(), mincvgArg.getValue(),
                       nthreadArg.getValue());

        }
        if (strcmp(argv[1], "summary") == 0){
            UnlabeledValueArg<string> sclustfileArg("sclustfile", "path of sclust file", true, "", "sclustfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            ValueArg<double> minlogLRArg("l","minlogLR","minimal logLR between joint and marigional probability, default: 0",
                                         false , 0, "minlogLR", cmd);
            ValueArg<int> mincountArg("c","mincount","minimal count of variants: 10", false , 10, "mincount", cmd);
            ValueArg<int> mincvgArg("v","cvg","minimal coverage of subspace: 10", false , 10, "mincvg", cmd);
            
            
            cmd.parse(argv2);
            
            SClust sclust;
            sclust.summary(sclustfileArg.getValue(), outfileArg.getValue(), minlogLRArg.getValue(),
                           mincountArg.getValue(), mincvgArg.getValue());
        }
        
        if (strcmp(argv[1], "eval") == 0){
            UnlabeledValueArg<string> patternfileArg("patternfile", "path of pattern file", true, "", "patternfile", cmd);
            UnlabeledValueArg<string> snpfileArg("snpfile", "path of snp file", true, "", "snpfile", cmd);
            
            cmd.parse(argv2);
            
            string out_file = patternfileArg.getValue() + ".eval";
            SClust sclust;
            sclust.eval_pattern(patternfileArg.getValue(), snpfileArg.getValue(), out_file);
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
            
            //SwitchArg isfastArg("q", "fast", "use fast algorithm to run dforest", cmd, false);
            //SwitchArg isinterArg("i", "intermediate", "output intermediate results", cmd, false);
            
            cmd.parse(argv2);
            cout << "minreads = " << minreadsArg.getValue() << endl;
            cout << "maxdepth = " << maxdepthArg.getValue() << endl;
            cout << "minfreq = " << minfreqArg.getValue() << endl;
            cout << "nthread = " << nthreadArg.getValue() << endl;
            
            AlignReaderM5 alignreader;
            AlignCoderSNV aligncoder;
            DForestSNVMax forestsnvmax(&alignreader, &aligncoder);
            //DForestSNV forestsnv(&alignreader, &aligncoder);
            //DForestSNVFast forestsnvfast(&alignreader, &aligncoder);
            
            DForest *ptr_forest;
            ptr_forest = &forestsnvmax;
            /*if (isfastArg.getValue()){
                ptr_forest = &forestsnvfast;
            }else{ 
                if (!isinterArg.getValue()){
                    ptr_forest = &forestsnvmax;
                }else{
                    ptr_forest = &forestsnv;
                }
            }*/
            string shell_cmd = "mkdir -p " + tmpdirArg.getValue();
            cout << shell_cmd << endl;
            system(shell_cmd.c_str());
            
            ptr_forest->run(encodefileArg.getValue(), alignfileArg.getValue(), cmpreadsfileArg.getValue(), outfileArg.getValue(), tmpdirArg.getValue(), minreadsArg.getValue(), maxdepthArg.getValue(), nthreadArg.getValue(), minfreqArg.getValue());
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
            ValueArg<double> mincondprobArg("c","condprob","minimal conditional probability, default: 0.8", false , 0.8, "condprob", cmd);
            
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
            
            cmd.parse(argv2);
            ErrorModelSNV errormodel;
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

        
        // mask encode file
        if (strcmp(argv[1], "mask")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> regionfileArg("regionfile", "path of region file", true, "", "regionfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            SwitchArg is0basedArg("b", "0based", "is 0-based", cmd, false);
            
            cmd.parse(argv2);
            
            AlignCoderSNV aligncodersnv;
            HClust hclust(&aligncodersnv);
            hclust.mask(encodefileArg.getValue(), regionfileArg.getValue(), outfileArg.getValue(), is0basedArg.getValue());
            
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

        if (strcmp(argv[1], "check_contained_reads")==0) {
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            //UnlabeledValueArg<string> outfileArg("outfile", "path of output files", true, "", "outfile", cmd);
            
            cmd.parse(argv2);
            string encode_file = encodefileArg.getValue();
            string align_file = alignfileArg.getValue();
            
            cout << "load encode_data and reads_range" << endl;
            vector<vector<int> > encode_data; loadencodedata(encode_data, encode_file);
            vector<ReadRange> reads_range; loadreadsrange(reads_range, align_file);
            
            cout << "check non-contained reads" << endl;
            Assembler assembler;
            vector<int> read_sel_idx = assembler.check_contained_reads(encode_data, reads_range);
            
            cout << "output results" << endl;
            select_lines(read_sel_idx, encode_file, encode_file + ".non_contained");
            select_lines(read_sel_idx, align_file, align_file + ".non_contained");

            string outfile = encode_file + ".non_contained.idx";
            ofstream fs_outfile;
            open_outfile(fs_outfile, outfile);
            for (int i=0; i < (int) read_sel_idx.size(); ++i)
                fs_outfile << read_sel_idx[i] << endl;
            fs_outfile.close();
            
            outfile = encode_file + ".non_contained.check_follower";
            vector<int> idx_with_follower = assembler.find_follower_reads(encode_data, reads_range, read_sel_idx, outfile);
            select_lines(idx_with_follower, encode_file, encode_file + ".non_contained.with_follower");
            select_lines(idx_with_follower, align_file, align_file + ".non_contained.with_follower");
        }
        
        if (strcmp(argv[1], "olc") == 0){
            // parse arguments
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            ValueArg<int> candsizeArg("s","candsize","maximal candidate size, default: 5", false , 5, "candsize", cmd);
            ValueArg<int> mincountArg("c","mincount","minimal count of variants: 10", false , 10, "mincount", cmd);
            ValueArg<double> mincondprobrg("p","minCondProb","minimal conditional probability, default: 0.15", false , 0.15, "minCondProb", cmd);
            ValueArg<double> maxcondprobrg("q","maxCondProb","maximal conditional probability, default: 0.75", false , 0.75, "maxCondProb", cmd);
            
            ValueArg<int> minmatchArg("n","minmatch","minimal number of matches between two reads, default: 2", false , 1, "minmatch", cmd);
            ValueArg<double> minsimArg("m","minsim","minimal similarity between two reads, default: 0.7", false , 0.7, "minsim", cmd);
            ValueArg<double> minpropArg("k","minprop","minimal proportion of match, default: 0.5", false , 0.5, "minprop", cmd);
            SwitchArg iscontainArg("r", "contain", "is check contained reads, default : false", cmd, false);
            
            cmd.parse(argv2);
            Assembler assembler(candsizeArg.getValue(), 20, mincountArg.getValue(),
                                mincondprobrg.getValue(), maxcondprobrg.getValue());
            assembler.olc(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(), minmatchArg.getValue(), minsimArg.getValue(), minpropArg.getValue(), iscontainArg.getValue());
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

        // correct reads
        if (strcmp(argv[1], "correct") == 0){
            // parse arguments
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> cmpreadsfileArg("cmpreads_diff_file", "path of cmpreads_diff file", true, "", "cmpreads_diff_file", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            
            
            ValueArg<int> candsizeArg("s","candsize","maximal candidate size, default: 5", false , 5, "candsize", cmd);
            ValueArg<int> mincountArg("c","mincount","minimal count of variants: 10", false , 10, "mincount", cmd);
            ValueArg<double> mincondprobrg("p","minCondProb","minimal conditional probability, default: 0.15", false , 0.15, "minCondProb", cmd);
            ValueArg<double> maxcondprobrg("q","maxCondProb","maximal conditional probability, default: 0.75", false , 0.75, "maxCondProb", cmd);

            cmd.parse(argv2);
            Assembler assembler(candsizeArg.getValue(), 20, mincountArg.getValue(),
                                mincondprobrg.getValue(), maxcondprobrg.getValue());
            vector<int> read_ids = assembler.correct_reads(encodefileArg.getValue(), alignfileArg.getValue(), cmpreadsfileArg.getValue(), outfileArg.getValue(), true);
            select_lines(read_ids, alignfileArg.getValue(), outfileArg.getValue() + ".m5");
            
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
