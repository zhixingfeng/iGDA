//
//  test_sclust.cpp
//  iGDA
//
//  Created by Zhixing Feng on 7/17/17.
//  Copyright (c) 2017 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/sclust/sclust.h"

#include <ctime>

bool cmpreads_topn_test(string encode_file, string align_file, string out_file, int topn = 10, double min_overlap = 0.25,
                          bool is_rm_single=true, bool is_binary=true, bool is_print_read_id=false, bool is_condprob=true)
{
    // load encode data
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);
    
    // load reads range
    vector<ReadRange> reads_range;
    loadreadsrange(reads_range, align_file);
    
    if (encode_data.size() != reads_range.size())
        throw runtime_error("cmpreads: size of encode_data and reads_range do not match.");
    
    cout << encode_data.size() << endl;
    
    // get the right-most variant location to determing size of template array
    int temp_array_size = 0;
    for (int i = 0; i < encode_data.size(); i++)
        for (int j = 0; j < encode_data[i].size(); j++)
            if (encode_data[i][j] > temp_array_size)
                temp_array_size = encode_data[i][j];
    temp_array_size++;
    vector<int> temp_array(temp_array_size, -1);
    
    // open output file
    FILE *p_out_file = NULL;
    if (is_binary){
        p_out_file = fopen(out_file.c_str(), "wb");
        if (p_out_file==NULL)
            throw runtime_error("fail to open out_file");
    }
    else{
        p_out_file = fopen(out_file.c_str(), "w");
        if (p_out_file==NULL)
            throw runtime_error("fail to open out_file");
        
    }
    
    // pairwise comparison
    for (int i=0; i<(int)encode_data.size(); i++){
        if ((i+1)%1000==0) cout << i+1 << endl;
        
        // fill the template array by the variants in the ith read
        for (int j = 0; j < encode_data[i].size(); j++)
            temp_array[encode_data[i][j]] = i;
        
        // store matches of the jth read to the ith read
        vector<ReadMatch> the_matches (encode_data.size(), ReadMatch());
        
        // compare other reads to cur_variant
        for (int j=0; j<(int)encode_data.size(); j++){
            if (j == i)
                continue;
            
            // get size of overlap of the two reads
            the_matches[j].start = reads_range[i].first > reads_range[j].first ? reads_range[i].first : reads_range[j].first;
            the_matches[j].end = reads_range[i].second < reads_range[j].second ? reads_range[i].second : reads_range[j].second;
            int n_overlap = the_matches[j].end - the_matches[j].start + 1;
            
            if (n_overlap < min_overlap * (reads_range[i].second - reads_range[i].first + 1))
                continue;
            
            
            // get intersection between two reads
            vector<int> cur_match;
            for (int k = 0; k < encode_data[j].size(); k++)
                if (temp_array[encode_data[j][k]] == i)
                    cur_match.push_back(encode_data[j][k]);
            
            
            the_matches[j].matches = cur_match;
            the_matches[j].n_overlap = n_overlap;
            
            if (is_condprob){
                if (encode_data[j].size() > 0)
                    the_matches[j].match_rate = (double) cur_match.size() / encode_data[i].size();
                else
                    the_matches[j].match_rate = 0;
            }else{
                the_matches[j].match_rate = (double) cur_match.size() / n_overlap;
            }
        }
        
        // sort the_matches according to match_rate
        stable_sort(the_matches.begin(), the_matches.end(), [](const ReadMatch & dl, const ReadMatch & dr) {return dl.match_rate > dr.match_rate;});
        
        // print topn matches
        //cout << "topn = " << topn << endl;
        //cout << "the_matches.size() = " << the_matches.size() << endl;
        int cur_size = topn <= the_matches.size() ? topn : (int)the_matches.size();
        //cout << "cur_size = " << cur_size << endl;
        for (int j=0; j<cur_size; j++){
            // skip matches with size < 2 if is_rm_single is true
            if (is_rm_single){
                if (the_matches[j].matches.size() < 2)
                    continue;
            }else{
                if (the_matches[j].matches.size() == 0)
                    continue;
            }
            // print results
            if (is_binary){
                int cur_match_size = (int)the_matches[j].matches.size();
                if (is_print_read_id){
                    fwrite(&i, sizeof(int), 1, p_out_file);
                    fwrite(&the_matches[j].start, sizeof(int), 1, p_out_file);
                    fwrite(&the_matches[j].end, sizeof(int), 1, p_out_file);
                }
                fwrite(&cur_match_size, sizeof(int), 1, p_out_file);
                fwrite(&the_matches[j].matches[0], sizeof(int), cur_match_size, p_out_file);
            }else{
                if (is_print_read_id){
                    fprintf(p_out_file, "%d\t", i);
                    fprintf(p_out_file, "%d\t", the_matches[j].start);
                    fprintf(p_out_file, "%d\t", the_matches[j].end);
                }
                for (int k=0; k<(int)the_matches[j].matches.size(); k++)
                    fprintf(p_out_file, "%d,", the_matches[j].matches[k]);
                fprintf(p_out_file, "\n");
            }
            
            
        }
        
    }
    cout << encode_data.size() << endl;
    
    fclose(p_out_file);
    
    return true;
}


double test_cal_nlogn(double x)
{
    if (int(x) == 0)
        return 0;
    return x*log(x);
}

double test_cal_logL_H1(double n_11, double n_10, double n_01, double n_00)
{
    double N = n_11 + n_10 + n_01 + n_00;
    double p_11 = n_11 / N; double p_10 = n_10 / N; double p_01 = n_01 / N; double p_00 = n_00 / N;
    double logLR_1 = n_11*log(p_11) + n_10*log(p_10) + n_01*log(p_01) + n_00*log(p_00);

    return logLR_1;
}

double test_cal_logL_H0(double n_11, double n_10, double n_01, double n_00)
{
    double N = n_11 + n_10 + n_01 + n_00;

    double p_11 = (n_11 + n_10)*(n_11 + n_01)/(N*N);
    double p_10 = (n_11 + n_10)/N - p_11;
    double p_01 = (n_11 + n_01)/N - p_11;
    double p_00 = 1 - p_11 - p_10 - p_01;
    double logLR_0 = n_11*log(p_11) + n_10*log(p_10) + n_01*log(p_01) + n_00*log(p_00);
    
    return logLR_0;
}


double test_cal_logLR(double n_11, double n_10, double n_01, double n_00){
    // log LR
    return test_cal_logL_H1(n_11, n_10, n_01, n_00) - test_cal_logL_H0(n_11, n_10, n_01, n_00);
}

TEST_CASE("test sclust::run(), debug", "[hide]")
{
    string align_file = "../results/sclust/align/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.clean.m5";
    string encode_file = "../results/sclust/encode/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.encode";
    string cmpreads_file = "../results/sclust/cmpreads/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.condprob.s15.cmpreads.20";
    
    string out_file = "../results/sclust/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n1.readid.1_vs_rest.20.debug";
    
    SClust sclust;
    sclust.run(encode_file, align_file, cmpreads_file, out_file, "./", 15, 0, 0, 0, 1);
}


TEST_CASE("test sclust::run()", "[hide]")
{
    string align_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string cmpreads_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.cmpreads.split";
    
    string out_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.sclust";
    
    SClust sclust;
    //sclust.run(encode_file, align_file, cmpreads_file, out_file, "./", 15, 0, 0, 0, 1);
    sclust.run(encode_file, align_file, cmpreads_file, out_file, "./", 5, 0, 10, 10, 1);
}

TEST_CASE("test sclust::run(), multithread", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads_topn_readid.bin";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n4.sclust";
    
    SClust sclust;
    sclust.run(encode_file, align_file, cmpreads_file, out_file, "../results/tmp", 15, 0, 0, 0, 4);
}

TEST_CASE("test sclust::eval_pattern", "[hide]")
{
    string pattern_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n8.summary.pattern.uniq.100";
    string true_snp_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.true.snp";
    string out_file = pattern_file + ".eval";
    SClust sclust;
    sclust.eval_pattern(pattern_file, true_snp_file, out_file);
}


TEST_CASE("test sclust::cal_logLR()","[hide]")
{
    SClust sclust;
    int n_11 = 100; int n_10 = 194; int n_01 = 493; int n_00 = 3982;
    double logLR = sclust.cal_logLR (n_11, n_10, n_01, n_11 + n_10 + n_01 + n_00);
    double test_logLR = test_cal_logLR (n_11, n_10, n_01, n_00);
    cout << "log LR: " << logLR << endl;
    cout << "test log LR: " << test_logLR << endl;
}

TEST_CASE("Test sclust::summary()", "[hide]")
{
    //string sclust_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n1.readid.10000";
    //string out_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n1.readid.10000.summary";
    //string sclust_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n1.readid.1_vs_rest.10000";
    //string out_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n1.readid.1_vs_rest.10000.summary";
    string sclust_file = "../results/B_10_cons_out_topn_dforestmax_n1.1_vs_rest.range.sclust";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n1.1_vs_rest.range.sclust.summary.subspace";
    SClust sclust;
    sclust.summary(sclust_file, out_file, 20, 10, 10);
    //sclust.summary(pattern_file, out_file, 1, 20);

}



TEST_CASE("Test sclust::split_subspace()", "[hide]")
{
    string cmpreads_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.cmpreads";
    string out_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.cmpreads.split";
    SClust sclust;
    sclust.split_subspace(cmpreads_file, out_file, 5);
    
}

TEST_CASE("Test cmpreads_topn, first 5000")
{
    string align_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000";
    string encode_file = "../data/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000";
    string cmpreads_file = "../results/detect_comb/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.5000.cmpreads";
    
    cmpreads_topn_test(encode_file, align_file, cmpreads_file, 20, 0.5,
                  true, true, true, true);
}






