//
//  sort_order.h
//  iGDA
//
//  Created by Zhixing Feng on 16/12/5.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef sort_order_h
#define sort_order_h

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class COMP_INT_GREATER
{
public:
    COMP_INT_GREATER(vector<int> &x) : base_arr(&x) {}
    virtual ~COMP_INT_GREATER(){}
    bool operator()(int a, int b)
    {
        return ((*base_arr)[a] > (*base_arr)[b]);
    }
private:
    vector <int> *base_arr;
};

class COMP_INT_LESS
{
public:
    COMP_INT_LESS(vector<int> &x) : base_arr(&x) {}
    virtual ~COMP_INT_LESS(){}
    bool operator()(int a, int b)
    {
        return ((*base_arr)[a] < (*base_arr)[b]);
    }
private:
    vector <int> *base_arr;
};

class COMP_DOUBLE_GREATER
{
public:
    COMP_DOUBLE_GREATER(vector<double> &x) : base_arr(&x) {}
    virtual ~COMP_DOUBLE_GREATER(){}
    bool operator()(double a, double b)
    {
        return ((*base_arr)[a] > (*base_arr)[b]);
    }
private:
    vector <double> *base_arr;
};

class COMP_DOUBLE_LESS
{
public:
    COMP_DOUBLE_LESS(vector<double> &x) : base_arr(&x) {}
    virtual ~COMP_DOUBLE_LESS(){}
    bool operator()(double a, double b)
    {
        return ((*base_arr)[a] < (*base_arr)[b]);
    }
private:
    vector <double> *base_arr;
};



inline vector<int> sort_order(vector<int> &x, bool is_decreasing=false)
{
    COMP_INT_GREATER comp_int_greater(x);
    COMP_INT_LESS comp_int_less(x);

    vector<int> idx(x.size(), 0);
    for (int i = 0; i < idx.size(); i++)
        idx[i] = i;
    if (is_decreasing)
        sort(idx.begin(), idx.end(), comp_int_greater);
    else
        sort(idx.begin(), idx.end(), comp_int_less);
    return idx;
}

inline vector<int> sort_order(vector<double> &x, bool is_decreasing=false)
{
    COMP_DOUBLE_GREATER comp_double_greater(x);
    COMP_DOUBLE_LESS comp_double_less(x);
    
    vector<int> idx(x.size(), 0);
    for (int i = 0; i < idx.size(); i++)
        idx[i] = i;
    if (is_decreasing)
        sort(idx.begin(), idx.end(), comp_double_greater);
    else
        sort(idx.begin(), idx.end(), comp_double_less);
    return idx;
}

#endif /* sort_order_h */
