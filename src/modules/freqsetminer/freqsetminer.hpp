//
//  freqsetminer.hpp
//  iGDA
//
//  Created by Zhixing Feng on 16/8/18.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef freqsetminer_hpp
#define freqsetminer_hpp

#include <headers.h>

class FreqSetMiner
{
public:
    FreqSetMiner(){}
    virtual ~FreqSetMiner(){}
    
protected:
    bool readline();
};



#endif /* freqsetminer_hpp */
