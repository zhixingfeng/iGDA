//
//  test_performance.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/11/22.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/misc.h"
TEST_CASE("compare speed of unordered_map vs unordered_set vs direct array search", "[hide]"){
    // conclusion: speed of array construction is >10x faster than hash; speed of array access is >60x faster than hash
    int N = 1;
    int n = 1048576;
    
    // speed of array allocation
    clock_t t_begin = clock();
    int t = 0;
    for (long int i=0; i<N; i++){
        vector<int> x(n,1);
        t += x[i];
    }
    cout << t << endl;
    clock_t t_end = clock();
    cout << "time for array allocation : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of unordered_map construction
    t_begin = clock();
    for (int i=0; i<N; i++){
        unordered_map<int, int> x;
        for (int j=0; j<n; j++)
            x[j] = 0;
    }
    t_end = clock();
    cout << "time for hash table construction : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of unordered_set construction
    vector<int> tmp(n,1);
    t_begin = clock();
    t = 0;
    for (int i=0; i<N; i++){
        unordered_set <int> x(tmp.begin(), tmp.end());
    }
    t_end = clock();
    cout << "time for unordered_set construction : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;

    
    // speed of array access
    int *x = new int[n];
    t_begin = clock();
    
    for (int i=0; i<N; i++){
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
    for (int i=0; i<N; i++){
        for (int j=0; j<n; j++)
            int y = z[j];
    }
    t_end = clock();
    cout << "time for unordered_map access : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // speed of unordered_set access
    unordered_set<int> q(tmp.begin(), tmp.end());
    unordered_set<int>::iterator y;
    t_begin = clock();
    for (int i=0; i<N; i++){
        for (int j=0; j<n; j++){
            y = q.find(j);
        }
    }
    if (y!=q.end())
        cout << *y << endl;
    
    t_end = clock();
    cout << "time for unordered_set access : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}

TEST_CASE("compare array and unordered_set iteration","[hide]")
{
    int n = 100000;

    // array iteration
    vector<int> x_vec(n,0);
    int y;
    clock_t t_begin = clock();
    for (int i=0; i<n; i++){
        x_vec[i] = i;
    }
    clock_t t_end = clock();
    cout << "time for array iteration : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;

    // unordered_set construction
    t_begin = clock();
    unordered_set<int> x_set;
    for (int i=1; i < n; i++)
        x_set.insert(i);
    t_end = clock();
    cout << "time for unordered_set construction : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // unordered_set iteration
    t_begin = clock();
    unordered_set<int>::iterator it;
    for (it=x_set.begin(); it!=x_set.end(); it++)
        y = *it;
    t_end = clock();
    cout << "time for unordered_set iteration : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;

}

TEST_CASE("sparse vetor iteration", "[hide]"){
    int n = 1000000;
    // array iteration
    vector<int> x_vec(n,0);
    x_vec[n/4] = 1;
    x_vec[2*n/4] = 2;
    x_vec[3*n/4] = 3;
    x_vec[n-1] = 4;
    
    vector<int> y_vec;
    clock_t t_begin = clock();
    for (int i=0; i<n; i++){
        if (x_vec[i]!=0)
            x_vec[i]++;
        //y_vec.push_back(x_vec[i]);
    }
    clock_t t_end = clock();
    cout << "time for array iteration : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
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

TEST_CASE("test file reading speed","[hide]")
{
    clock_t t_begin = clock();
    std::ifstream in("../results/MSSA_61_forward_encode_snv_cmpreads_col2.txt", std::ios::in | std::ios::binary);
    std::string contents;
    if (in)
    {
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
    }
    clock_t t_end = clock();
    cout << "size of file: " << contents.size() << endl;
    cout << "time of read(): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    t_begin = clock();
    in.open("../results/MSSA_61_forward_encode_snv_cmpreads_col2.txt", std::ios::in | std::ios::binary);
    while (!in.eof()){
        getline(in, contents);
        vector<int> cur_data = split_int(contents,',');
    }
    t_end = clock();
    cout << "time of getline(): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}

/*TEST_CASE("Test file writing speed", "[hide]"){
    // fprintf is > 5 times faster than "<<" on my MacPro (both true for -O3 and non-optimized version)
    DforestResult rl;
    string outfile_c = "../results/test_output_c.txt";
    string outfile_cpp = "../results/test_output_cpp.txt";
    
    int B = 1000000;
    // C style writing
    clock_t t_begin = clock();
    FILE *p_c = fopen(outfile_c.c_str(), "w");
    if (p_c == NULL)
        throw runtime_error("unable to open outfile_c");
    for (int i = 0; i < B; i++){
        fprintf(p_c, "%lf,%lf,%d,%d\n", rl.bf, rl.p_y_xp, rl.n_y_xp, rl.n_xp);
    }
    fclose(p_c);
    clock_t t_end = clock();
    cout << "time of c style writing: " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    
    // Cpp stype writing 
    t_begin = clock();
    ofstream p_cpp; open_outfile(p_cpp, outfile_cpp);
    for (int i = 0; i < B; i++)
        p_cpp<<rl.bf <<',' << rl.p_y_xp <<',' << rl.n_y_xp <<',' << rl.n_xp << endl;
    p_cpp.close();
    t_end = clock();
    cout << "time of cpp style writing: " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;

}*/
TEST_CASE("Test size of data type","[hide]")
{
    cout << "char : " << sizeof(char) << endl;
    cout << "int : " << sizeof(int) << endl;
    cout << "long : " << sizeof(long) << endl;
    cout << "long int : " << sizeof(long int) << endl;
    cout << "double : " << sizeof(double) << endl;
    
}


TEST_CASE("Test convert text to binary file","[hide]")
{
    ifstream in; 
    open_infile(in, "../results/MSSA_61_forward_encode_snv_cmpreads_col2.txt");
    FILE * p_outfile = fopen("../results/MSSA_61_forward_encode_snv_cmpreads_col2.bin", "wb");
    if (p_outfile==NULL)
        runtime_error("fail to open outfile");
    
    string buf;
    while(!in.eof()){
        getline(in, buf);
        if (in.eof())
            break;
        vector<int> x = split_int(buf, ',');
        int y = (int)x.size();
        fwrite(&y, sizeof(int), 1, p_outfile);
        fwrite(&x[0], sizeof(int), (int)x.size(), p_outfile);
    }
    in.close();
    fclose(p_outfile);
}


TEST_CASE("Test read binary file", "[hide]")
{
    clock_t t_begin = clock();
    FILE *p_infile = fopen("../results/B_10_cons_cmpreads.bin", "rb");
    if (p_infile==NULL)
        runtime_error("fail to open infile");
    ofstream p_outfile;
    open_outfile(p_outfile, "../results/B_10_cons_cmpreads.bin.rev.txt");
    while(1){
        int cur_size ;
    
        fread(&cur_size, sizeof(int), 1, p_infile);
        vector<int> cur_data_vec(cur_size);
        fread(&cur_data_vec[0], sizeof(int), cur_size, p_infile);
        if (feof(p_infile))
            break;
        for (int i=0; i<cur_data_vec.size(); i++)
            p_outfile << cur_data_vec[i] << ',';
        p_outfile << endl;
    }
    
    fclose(p_infile);
    p_outfile.close();
    clock_t t_end = clock();
    cout << "time of read binary(): " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;

}


TEST_CASE("bit shift vs multiply","[hide]")
{
    int N = 1000000000;
    
    clock_t t_begin = clock();
    int x = 10;
    int y = 0;
    for (int i=0; i<N; i++){
        y = x << 6;
    }
    clock_t t_end = clock();
    cout << "time of bit shift: " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    cout << "y bit shift: " << y << endl;
    
    t_begin = clock();
    x = 10;
    for (int i=0; i<N; i++){
        y = x * 64;
    }
    t_end = clock();
    cout << "time of multiply: " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
    cout << "y multiply: " << y << endl;
}


TEST_CASE("test maximal size of vector", "[hide]")
{
    vector<int> myvector;
    cout << "maximal size of vector:" <<  myvector.max_size() << endl;
}




