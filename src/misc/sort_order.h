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

static vector <int> *base_arr_int;
static vector <double> *base_arr_double;


inline bool compar_less_int (int a, int b)
{
    return ((*base_arr_int)[a] < (*base_arr_int)[b]);
}

inline bool compar_greater_int (int a, int b)
{
    return ((*base_arr_int)[a] > (*base_arr_int)[b]);
}

inline bool compar_less_double (int a, int b)
{
    return ((*base_arr_double)[a] < (*base_arr_double)[b]);
}

inline bool compar_greater_double (int a, int b)
{
    return ((*base_arr_double)[a] > (*base_arr_double)[b]);
}

inline vector<int> sort_order(vector<int> &x, bool is_decreasing=false)
{
    base_arr_int = &x;
    vector<int> idx(x.size(), 0);
    for (int i = 0; i < idx.size(); i++)
        idx[i] = i;
    if (is_decreasing)
        sort(idx.begin(), idx.end(), compar_greater_int);
    else
        sort(idx.begin(), idx.end(), compar_less_int);
    return idx;
}

inline vector<int> sort_order(vector<double> &x, bool is_decreasing=false)
{
    base_arr_double = &x;
    vector<int> idx(x.size(), 0);
    for (int i = 0; i < idx.size(); i++)
        idx[i] = i;
    if (is_decreasing)
        sort(idx.begin(), idx.end(), compar_greater_double);
    else
        sort(idx.begin(), idx.end(), compar_less_double);
    return idx;
}

#endif /* sort_order_h */
