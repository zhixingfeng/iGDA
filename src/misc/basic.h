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



#endif /* basic_h */
