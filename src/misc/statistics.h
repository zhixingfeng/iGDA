//
//  statistics.h
//  iGDA
//
//  Created by Zhixing Feng on 10/4/16.
//  Copyright (c) 2016 Zhixing Feng. All rights reserved.
//

#ifndef iGDA_statistics_h
#define iGDA_statistics_h

#include <cmath>
#include "../tools/prob/prob.hpp"
#include <boost/math/distributions/beta.hpp>
#include "../../include/headers.h"
inline double lbeta(double a, double b)
{
    return boost::math::lgamma(a) + boost::math::lgamma(b) - boost::math::lgamma(a + b) ;
}


inline double pnorm(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}

// approximation of log Gamma function according to "Gergo Nemes, New asymptotic expansion for the Gamma function, Arch. Math. 95 (2010), 161–169"
inline double lgamma(double x)
{
    if (x<=2)
        return 0;
    else
        return (log(2*GDA_PI) - log(x))/2 + x*( log( x + 1/(12*x-1/(10*x)) ) - 1 );
}

// binomial log likelihood ratio
inline double binom_log_lr (double x, double n, double p0)
{
    double p1 = x/n;
    double log_lr = x*( log(p1) - log(p0) ) + (n-x)*( log(1-p1) - log(1-p0) );
    return log_lr;
}

// binomial log bayes factor
inline double binom_log_bf(double x, double n, double p0)
{
    double a = 1;
    double b = 1;
    
    //double H1_part1 = r8_beta(x+a, n-x+b);
    double H1_part1 = boost::math::beta(x+a, n-x+b);
    H1_part1 = H1_part1 < EPS ? EPS : H1_part1;
    
    //double H1_part2 = 1 - beta_cdf(p0, x+a, n-x+b);
    double H1_part2 = boost::math::ibetac(x+a, n-x+b, p0);
    H1_part2 = H1_part2 < EPS ? EPS : H1_part2;
    
    double log_H1 = log(H1_part1) + log(H1_part2);
    //double log_H1 = log(r8_beta(x+a, n-x+b)) + log(1 - beta_cdf(p0, x+a, n-x+b));
    double log_H0 = x*log(p0) + (n-x)*log(1-p0);
    double log_bf = log_H1 - log_H0;
    
    return log_bf;
}

// binomial log bayes factor (legacy, actually it is a two-side test not an exact one-side test)
inline double binom_log_bf_legacy (double x, double n, double p0)
{
    double L_1 = lgamma(x+1) + lgamma(n-x+1) - lgamma(n+2);
    double L_0 = x*log(p0) +(n-x)*log(1-p0);
    
    double log_bf = L_1 - L_0;
    // restrict minimal bf to be 0 and x<n*p0 make sure it is one-side test;
    if (log_bf<0 || n==0 || x<n*p0)
        log_bf=0;
    return log_bf;
}

#endif
