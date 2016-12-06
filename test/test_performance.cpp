//
//  test_performance.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/11/22.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"


TEST_CASE("compare speed of unordered_map vs unordered_set vs direct array search","[hide]"){
    // conclustion: speed of array construction is >10x faster than hash; speed of array access is >60x faster than hash
    int N = 50000;
    int n = 10000;
    
    // speed of array allocation
    clock_t t_begin = clock();
    int t = 0;
    for (long int i=0; i<10000; i++){
        vector<int> x(N,1);
        t += x[i];
    }
    cout << t << endl;
    clock_t t_end = clock();
    cout << "time for array allocation : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of unordered_map construction
    t_begin = clock();
    for (int i=0; i<10000; i++){
        unordered_map<int, int> x;
        for (int j=0; j<n; j++)
            x[j] = 0;
    }
    t_end = clock();
    cout << "time for hash table construction : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of unordered_map construction
    vector<int> tmp(n,1);
    t_begin = clock();
    t = 0;
    for (int i=0; i<10000; i++){
        unordered_set <int> x(tmp.begin(), tmp.end());
    }
    t_end = clock();
    cout << "time for unordered_set construction : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;

    
    // speed of array access
    int *x = new int[N];
    t_begin = clock();
    
    for (int i=0; i<10000; i++){
        for (int j=0; j<n; j++){
            int y = x[j];
            x[j] = 1;
        }
    }
    t_end = clock();
    x[0] += 1;
    delete[] x;
    
    cout << "time for array access : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of unordered_map access
    unordered_map<int, int> z;
    for (int j=0; j<n; j++)
        z[j] = 0;
    t_begin = clock();
    for (int i=0; i<10000; i++){
        for (int j=0; j<n; j++)
            int y = z[j];
    }
    t_end = clock();
    cout << "time for unordered_map access : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of unordered_set access
    unordered_set<int> q(tmp.begin(), tmp.end());
    unordered_set<int>::iterator y;
    t_begin = clock();
    for (int i=0; i<10000; i++){
        for (int j=0; j<n; j++){
            y = q.find(j);
        }
    }
    if (y!=q.end())
        cout << *y << endl;
    
    t_end = clock();
    cout << "time for unordered_set access : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}

TEST_CASE("Compare iteration speed of array and list", "[hide]"){
    int B = 1000;
    int N = 10000;
    
    struct Node {
        Node (): value(0),ptr(NULL){}
        int value;
        Node *ptr;
    };

    // create array
    
    vector<Node> x_arr(N,Node());
    //int x_arr[N];
    
    // access array
    clock_t t_begin = clock();
    for (long int i = 0; i < B; ++i){
        for (int j = 0; j < N; ++j){
            int x = x_arr[j].value;
            x_arr[j].value = 1;
            //int x += x_arr[j];
            //x_arr[j] = 1;
        }
    }
    clock_t t_end = clock();
    x_arr[0].value += 1;
    //x_arr[0] += 1;
    cout << "time to access array : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // create list
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
            int x = it->value;
            it->value = 1;
            it = it->ptr;
        }
    }
    t_end = clock();
    cout << "time to access list: " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
}

bool myfunction (int i,int j) { return (i<j); }

TEST_CASE("test speed of sort of std library", "[hide]")
{
    int B = 100000;
    std::vector<int> myvector;
    
    // set some values:
    for (int i=1; i<100; i++) myvector.push_back(i);   // 1 2 3 4 5 6 7 8 9

    mt19937 r{std::random_device{}()};
    shuffle(std::begin(myvector), std::end(myvector), r);

    clock_t t_begin = clock();
    for (int i=0; i<B; i++){
        shuffle(std::begin(myvector), std::end(myvector), r);
    }
    clock_t t_end = clock();
    clock_t t_shuffle = t_end - t_begin;
    cout << "time to shuffle vector: " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    t_begin = clock();
    for (int i=0; i<B; i++){
        sort(myvector.begin(), myvector.end());
        shuffle(std::begin(myvector), std::end(myvector), r);
    }
    t_end = clock();
    cout << "time to sort: " << double(t_end - t_begin - t_shuffle)/CLOCKS_PER_SEC << endl;
    
}

