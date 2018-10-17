//
//  test_prob.cpp
//  iGDA
//
//  Created by Zhixing Feng on 7/5/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../tools/prob/prob.hpp"
#include <boost/math/distributions/beta.hpp>

TEST_CASE("test cdf of beta distribution", "[hide]")
{
    cout << "beta_cdf(0.1, 2, 5) = " << beta_cdf(0.1, 2, 5) << endl;
}

TEST_CASE("test beta function", "[hide]")
{
    cout << "r8_beta(2, 5) = " << r8_beta(2, 5) << endl;
    cout << "r8_beta(1, 1) = " << r8_beta(1, 1) << endl;
}

TEST_CASE("test beta distribution of boost", "[hide]")
{
    cout << "ibeta(2,5,0.02) = " << boost::math::ibeta(2,5,0.02) << endl;
    cout << "ibetac(2,5,0.02) = " << boost::math::ibetac(2,5,0.02) << endl;
    cout << "ibeta(7200, 280000 - 7200 ,0.02) = " << boost::math::ibeta(7200, 280000 - 7200 ,0.02) << endl;
    cout << "ibetac(7200, 280000 - 7200 ,0.02) = " << boost::math::ibetac(7200, 280000 - 7200 ,0.02) << endl;
    
}

