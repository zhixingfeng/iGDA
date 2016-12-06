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

#endif /* basic_h */
