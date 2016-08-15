//
//  seqopt.h
//  iGDA
//
//  Created by Zhixing Feng on 16/8/15.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#ifndef seqopt_h
#define seqopt_h

#include "../../include/headers.h"
// get reverse complementary DNA sequence
inline string getrevcomp(string dnaseq)
{
    string revseq = "";
    int len = (int) dnaseq.size();
    for (int i=len-1; i>=0; i--) {
        char curbase = dnaseq[i];
        switch(curbase){
            case 'A':
                curbase = 'T';
                break;
            case 'C':
                curbase = 'G';
                break;
            case 'G':
                curbase = 'C';
                break;
            case 'T':
                curbase = 'A';
        }
        revseq.push_back(curbase);
    }
    return revseq;
}


#endif /* seqopt_h */
