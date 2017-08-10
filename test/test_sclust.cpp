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

TEST_CASE("test sclust::run()","[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads_topn_readid.bin";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n1.sclust";
    
    SClust sclust;
    //sclust.run(encode_file, align_file, cmpreads_file, out_file, "./", 15, 0, 0, 0, 1);
    sclust.run(encode_file, align_file, cmpreads_file, out_file, "./", 15, 0, 10, 10, 1);
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

TEST_CASE("Test sclust::summary()")
{
    string sclust_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n1.readid.10000";
    string out_file = "../results/sclust/ERR1109332_ERR1246962_ERR1246953_ERR1599920_ERR1588648.top20.sclust.lr.n1.readid.10000.summary";
    SClust sclust;
    sclust.summary(sclust_file, out_file, 20, 10, 10);
    //sclust.summary(pattern_file, out_file, 1, 20);

}









