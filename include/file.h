/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   file.h
 * Author: zhixingfeng
 *
 * Created on November 11, 2015, 5:50 PM
 */

#ifndef FILE_H
#define FILE_H
#include "stl.h"

inline void open_infile(ifstream &fs_infile, string filename){
    fs_infile.open(filename.c_str());
    if (!fs_infile.is_open()){
        string err_msg = "Error: fail to open " + filename;
        throw runtime_error(err_msg);
    }
}

inline void open_outfile(ofstream &fs_outfile, string filename){
    fs_outfile.open(filename.c_str());
    if (!fs_outfile.is_open()){
        string err_msg = "Error: fail to open " + filename;
        throw runtime_error(err_msg);
    }
}

#endif /* FILE_H */

