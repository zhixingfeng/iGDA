//
//  test_performance.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/11/22.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"


TEST_CASE("compare speed of hash table vs direct array search"){
    // conclustion: speed of array construction is >10x faster than hash; speed of array access is >60x faster than hash
    int N = 50000;
    int n = 20;
    
    // speed of array allocation
    clock_t t_begin = clock();
    for (int i=0; i<10000; i++){
        int *x = new int[N];
        delete[] x;
    }
    clock_t t_end = clock();
    cout << "time for array allocation : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of hash table construction
    t_begin = clock();
    for (int i=0; i<10000; i++){
        unordered_map<int, int> x;
        for (int j=0; j<n; j++)
            x[j] = 0;
    }
    t_end = clock();
    cout << "time for hash table construction : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of array access
    int *x = new int[N];
    t_begin = clock();
    for (int i=0; i<10000; i++){
        for (int j=0; j<n; j++)
            int y = x[j];
    }
    t_end = clock();
    delete[] x;
    cout << "time for array access : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of hash table access
    unordered_map<int, int> z;
    for (int j=0; j<n; j++)
        z[j] = 0;
    t_begin = clock();
    for (int i=0; i<10000; i++){
        for (int j=0; j<n; j++)
            int y = z[j];
    }
    t_end = clock();
    cout << "time for hash table access : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}

TEST_CASE("Compare iteration speed of array and list"){
    int B = 10000;
    int N = 10000;
    
    // create array
    int x_arr[N];
    
    
    // access array
    clock_t t_begin = clock();
    for (int i = 0; i < B; i++){
        for (int j = 0; j < N; j++){
            int x = x_arr[i];
        }
    }
    clock_t t_end = clock();
    cout << "time to access array : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // create list
    struct Node {
        Node (): value(0),ptr(NULL){}
        int value;
        Node *ptr;
    };
    Node *head_node = new Node;
    Node *pre_node = head_node;
    for (int i = 0; i < N; i++){
        Node *cur_node = new Node;
        pre_node->ptr = cur_node;
        pre_node = cur_node;
    }

    // access list
    t_begin = clock();
    for (int i = 0; i < B; i++){
        Node *it = head_node;
        while(it->ptr != NULL){
            //int x = it->value;
            it = it->ptr;
        }
    }
    t_end = clock();
    cout << "time to access list: " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
}



