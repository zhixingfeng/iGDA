//
//  basic.h
//  iGDA
//
//  Created by Zhixing Feng on 16/12/5.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef basic_h
#define basic_h

#include "../../include/headers.h"

// product of a vector
inline double prod(const vector<double> &x)
{
    if (x.size() == 0)
        throw runtime_error("product of an tempty vector");
    double y = 1;
    for (auto i = 0; i < x.size(); ++i)
        y *= x[i];
    return y;
}

// define "<<" to output vector
template<typename T>
std::ostream& operator<<(std::ostream& s, std::vector<T> t) { 
    for (std::size_t i = 0; i < t.size(); i++) {
        s << t[i] << (i == t.size() - 1 ? "" : ",");
    }
    return s;
}

// count the number of 1s in binary representation of the integer
inline int bitcount(uint32_t u)
{
    uint32_t uCount;
    
    uCount = u - ((u >> 1) & 033333333333) - ((u >> 2) & 011111111111);
    return ((uCount + (uCount >> 3)) & 030707070707) % 63;
}

// slide window for vector
template<typename T>
inline void slide_win(const vector<T> &vec_in, vector<T> &vec_out, int start, int win_size)
{
    if (start >= (int)vec_in.size())
        throw runtime_error("start >= (int)vec_in.size()");
    int end = start + win_size - 1 < vec_in.size() ? start + win_size - 1 : (int)vec_in.size() - 1;
    for (int i = start; i <= end; ++i)
        vec_out.push_back(vec_in[i]);
}

// intersect of sets
inline vector<int> intersect(const vector<int> &vec_1, const vector<int> &vec_2)
{
    unordered_set<int> vec_1_set(vec_1.begin(), vec_1.end());
    vector<int> vec_inter;
    for (auto i = 0; i < vec_2.size(); ++i){
        if (vec_1_set.find(vec_2[i]) != vec_1_set.end()){
            vec_inter.push_back(vec_2[i]);
        }
    }
    return vec_inter;
}


#endif /* basic_h */
