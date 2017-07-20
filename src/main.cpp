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
    cout << "command = samtofa, bamtofa, m5tofa, encode, cmpreads, sclust, bin2txt, txt2bin, dforest, sort, filter, contexteffect, merge, mergeall, mask, dist" << endl;
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
            SwitchArg istextArg("t", "text", "output text file", cmd, false);
            SwitchArg isdupArg("d", "dup", "keep duplicated candidates", cmd, false);
            cmd.parse(argv2);
            
            if (!isdupArg.getValue()){
                string tmpoutfile = outfileArg.getValue() + ".tmp";
                if (topnArg.getValue() > 0)
                    cmpreads_topn(encodefileArg.getValue(), alignfileArg.getValue(), tmpoutfile,
                                  topnArg.getValue(), overlapArg.getValue(), true, false);
                else
                    cmpreads(encodefileArg.getValue(), alignfileArg.getValue(), tmpoutfile, overlapArg.getValue(), true, false);
                if (istextArg.getValue()){
                    string shell_cmd = "awk '!seen[$0]++' " +  tmpoutfile + " > " + outfileArg.getValue();
                    cout << shell_cmd << endl;
                    system(shell_cmd.c_str());
                }else{
                    string shell_cmd = "awk '!seen[$0]++' " +  tmpoutfile + " > " + outfileArg.getValue() + ".uniq";
                    cout << shell_cmd << endl;
                    system(shell_cmd.c_str());
                    
                    shell_cmd = "igda txt2bin " + outfileArg.getValue() + ".uniq " + outfileArg.getValue();
                    cout << shell_cmd << endl;
                    system(shell_cmd.c_str());
                    
                    shell_cmd = "rm " + outfileArg.getValue() + ".uniq";
                    cout << shell_cmd << endl;
                    system(shell_cmd.c_str());
                }
                
                string shell_cmd = "rm " + tmpoutfile;
                cout << shell_cmd << endl;
                system(shell_cmd.c_str());

            }else{
                if (topnArg.getValue() > 0)
                    cmpreads_topn(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                                  topnArg.getValue(), overlapArg.getValue(), true, !istextArg.getValue());
                else
                    cmpreads(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(),
                             overlapArg.getValue(), true, !istextArg.getValue());
            }
        }
        
        // subspace clustering
        if (strcmp(argv[1], "sclust") == 0){
            UnlabeledValueArg<string> encodefileArg("encodefile", "path of encode file", true, "", "encodefile", cmd);
            UnlabeledValueArg<string> alignfileArg("alignfile", "path of align file", true, "", "alignfile", cmd);
            UnlabeledValueArg<string> cmpreadsfileArg("cmpreadsfile", "path of cmpreads file", true, "", "cmpreadsfile", cmd);
            UnlabeledValueArg<string> outfileArg("outfile", "path of output file", true, "", "outfile", cmd);
            UnlabeledValueArg<string> tmpdirArg("tmpdir", "temporary directory", true, "", "tmpdir", cmd);
            
            ValueArg<int> maxsubdimArg("s","maxsubdim","maximal subspace dimension, default: 15", false , 15, "maxsubdim", cmd);
            ValueArg<double> minratioArg("r","mincondprob","minimal conditional probability, default: 0.2", false , 0.2, "mincondprob", cmd);
            ValueArg<int> mincountArg("c","mincount","minimal count of variants: 10", false , 10, "mincount", cmd);
            ValueArg<int> mincvgArg("v","cvg","minimal coverage of subspace: 20", false , 20, "mincvg", cmd);
            ValueArg<int> nthreadArg("n","nthread","number of threads, default: 1", false , 1, "nthread", cmd);
            
            cmd.parse(argv2);
            
            SClust sclust;
            sclust.run(encodefileArg.getValue(), alignfileArg.getValue(), cmpreadsfileArg.getValue(),
                       outfileArg.getValue(), tmpdirArg.getValue(), maxsubdimArg.getValue(),
                       minratioArg.getValue(), mincountArg.getValue(), mincvgArg.getValue(),
                       nthreadArg.getValue());

        }
        // convert binary cmpreadsfile to text
        if (strcmp(argv[1], "bin2txt") == 0) {
            UnlabeledValueArg<string> binfileArg("binfile", "binary cmpreads file", true, "", "binfile", cmd);
            UnlabeledValueArg<string> txtfileArg("txtfile", "text cmpreads file", true, "", "txtfile", cmd);
            cmd.parse(argv2);
            
            cmpreads_bin2txt(binfileArg.getValue(), txtfileArg.getValue());
        }
        
        // convert text cmpreadsfile to binary
        if (strcmp(argv[1], "txt2bin")==0) {
            UnlabeledValueArg<string> txtfileArg("txtfile", "text cmpreads file", true, "", "txtfile", cmd);
            UnlabeledValueArg<string> binfileArg("binfile", "binary cmpreads file", true, "", "binfile", cmd);
            
            cmd.parse(argv2);
            
            cmpreads_txt2bin(txtfileArg.getValue(), binfileArg.getValue());
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
            
            SwitchArg isfastArg("q", "fast", "use fast algorithm to run dforest", cmd, false);
            SwitchArg isinterArg("i", "intermediate", "output intermediate results", cmd, false);
            
            cmd.parse(argv2);
            
            AlignReaderM5 alignreader;
            AlignCoderSNV aligncoder;
            DForestSNVMax forestsnvmax(&alignreader, &aligncoder);
            DForestSNV forestsnv(&alignreader, &aligncoder);
            DForestSNVFast forestsnvfast(&alignreader, &aligncoder);
            
            DForest *ptr_forest;
            if (isfastArg.getValue()){
                ptr_forest = &forestsnvfast;
            }else{ 
                if (!isinterArg.getValue()){
                    ptr_forest = &forestsnvmax;
                }else{
                    ptr_forest = &forestsnv;
                }
            }
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
            SwitchArg isnmissArg("n", "nmiss", "is ouput number of mismatches", cmd, false);
            
            cmd.parse(argv2);
            
            HClust hclust;
            hclust.dist(encodefileArg.getValue(), alignfileArg.getValue(), outfileArg.getValue(), isnmissArg.getValue());
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
