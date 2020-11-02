//
//  test_dforest.cpp
//  iGDA
//
//  Created by Zhixing Feng on 16/12/2.
//  Copyright © 2016年 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/modules/alignreader/alignreaderm5.h"
#include "../src/modules/aligncoder/aligncodersnv.h"
#include "../src/modules/dforest/dforestsnvmax.h"
#include "../src/modules/dforest/dforestsnvstl.h"
#include <ctime>
#include <stxxl/vector>
/*TEST_CASE("test DForest::run() (input from memory/stxxl)", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads_topn.txt";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n1.txt.stxxl";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVMax forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    cout << "load encode_data" << endl;
    vector<vector<int> > encode_data;
    loadencodedata(encode_data, encode_file);

    cout << "load alignment" << endl;
    AlignReaderM5 AlignReaderM5_obj;
    vector<Align> align_data;
    AlignReaderM5_obj.read(align_file, align_data);

    cout << "cmpreads" << endl;
    stxxl_vector_type_int cmpreads_data;
    loadcmpreads(cmpreads_data, cmpreads_file);

    cout << "dforest" << endl;
    int start_time= (int)clock();
    ptr_forest->run(encode_data, align_data, cmpreads_data, 8, 5, 1);
    int stop_time= (int)clock();
    
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
    
    // write results (unordered) to outfile
    ofstream fs_outfile;  open_outfile(fs_outfile, out_file);
    unordered_map<int, DforestResult> result = ptr_forest->get_result();
    for (auto it = result.begin(); it!=result.end(); ++it){
        if (it->second.link_loci.size() > 0 && it->second.p_y_xp >= 0){
            fs_outfile << it->second.focal_locus << '\t' << it->second.bf << '\t'
                        << it->second.p_y_xp << '\t' << it->second.n_y_xp << '\t'
                        << it->second.n_xp << '\t' << it->second.link_loci.size() << '\t';
            for (int j = 0; j < it->second.link_loci.size(); j++)
                fs_outfile << it->second.link_loci[j] << ',';
            fs_outfile << endl;
        }
    }
    fs_outfile.close();
    
}*/


/*TEST_CASE("test DForest::run()", "[hide]")
{
    string align_file = "../data/B_10_cons.m5";
    string encode_file = "../results/B_10_cons.encode";
    string cmpreads_file = "../results/B_10_cons_cmpreads_topn.bin";
    string out_file = "../results/B_10_cons_out_topn_dforestmax_n1.txt";
    //string out_file = "../results/dummy_cmpreads_out.txt";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    //DForestSNV forestsnv(&alignreader, &aligncoder);
    //DForestSNVFast forestsnv(&alignreader, &aligncoder);
    DForestSNVMax forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results" , 8, 5, 1);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}

TEST_CASE("test DForest::filter()", "[hide]")
{
    string dforest_file = "../results/B_10_cons.dforest";
    string out_file = "../results/B_10_cons.dforest.f0.5.filter";
    
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVMax forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    ptr_forest->filter(dforest_file, out_file, 0.5);
}*/


/*TEST_CASE("test DForest::build_tree()", "[hide]")
{
    string encode_file = "../results/B_10_cons.encode";
    string align_file = "../data/B_10_cons.m5";
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNV forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    ptr_forest->call_pileup_var(encode_file);
    ptr_forest->call_pileup_reads(align_file);
    
    vector<int> cand_loci({142, 556, 628, 702, 818, 1001, 1003, 1050, 10000});
    vector<vector<Result> > rl(ptr_forest->get_pileup_var().size(), vector<Result>());
    vector<int> temp_vec_var(ptr_forest->get_n_reads(), -1);
    vector<int> temp_vec_var_lock(ptr_forest->get_n_reads(), -1);
    vector<int> temp_vec_read(ptr_forest->get_n_reads(), -1);
    vector<int> temp_vec_read_lock(ptr_forest->get_n_reads(), -1);

    clock_t t_begin = clock();
    for (int i=0; i<1000; i++)
        ptr_forest->build_tree(cand_loci, temp_vec_var, temp_vec_var_lock, temp_vec_read, temp_vec_read_lock, 8, 10);
    clock_t t_end = clock();
    cout << "time for build_tree : " << double(t_end - t_begin)/CLOCKS_PER_SEC << endl;
}
*/


/*TEST_CASE("test DForest::run() (hbv)", "[hide]")
{
    string encode_file = "../results/LoFreq_HBV/igda_result_large/realign.encode";
    string align_file = "../results/LoFreq_HBV/igda_result_large/realign.m5";
    string cmpreads_file = "../results/LoFreq_HBV/igda_result_large/realign.cmpreads";
    string ref_file = "../results/LoFreq_HBV/cat_wild_large.fasta";
    string out_file = "../results/LoFreq_HBV/igda_result_large/realign.dforest.debug.v2";
    
    
    AlignReaderM5 alignreader;
    AlignCoderSNV aligncoder;
    DForestSNVSTL forestsnv(&alignreader, &aligncoder);
    DForest *ptr_forest = &forestsnv;
    
    ptr_forest->load_homo_blocks(ref_file);
    int start_time= (int)clock();
    ptr_forest->run(encode_file, align_file, cmpreads_file, out_file, "../results/LoFreq_HBV/igda_result_large/tmp" , 12, 10000, 1);
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}*/


TEST_CASE("test stxxl vector  speed", "[hide]")
{
    typedef stxxl::VECTOR_GENERATOR<int>::result vector_stxxl;
    
    cout << "compare push_back speed for stxxl vector and stl vector" << endl;
    
    // stxxl vector
    int start_time= (int)clock();
    vector_stxxl vector_1;
    for (int i = 0; i < 1024 * 1024 * 1024; i++)
    {
        vector_1.push_back(i + 2);
    }
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
    
    // stl vector
    start_time= (int)clock();
    vector<int> vector_2;
    for (int i = 0; i < 1024 * 1024 * 1024; i++)
    {
        vector_2.push_back(i + 2);
    }
    stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;

    
    cout << "compare random access speed for stxxl vector and stl vector" << endl;
    
    start_time= (int)clock();
    for (int i = 0; i < 1024 * 1024 * 1024; i++)
    {
        vector_1[i] -= 1;
    }
    stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
    
    start_time= (int)clock();
    for (int i = 0; i < 1024 * 1024 * 1024; i++)
    {
        vector_2[i] -= 1;
    }
    stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
}


TEST_CASE("test stxxl vector of vector", "[hide]")
{
    //typedef stxxl::VECTOR_GENERATOR<stxxl::VECTOR_GENERATOR<int>::result>::result vector_stxxl;

    //typedef vector<stxxl::VECTOR_GENERATOR<int>::result> vector_stxxl;
    typedef stxxl::vector<vector<int> > vector_stxxl;
    // stxxl vector
    int start_time= (int)clock();
    vector_stxxl vector_1;
    for (int i = 0; i < 1024 * 1024 ; i++)
    {
        vector_1.push_back(vector<int>(8,4));
    }
    int stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
    
    // stl vector
    start_time= (int)clock();
    vector<vector<int> > vector_2;
    for (int i = 0; i < 1024 * 1024 ; i++)
    {
        vector_2.push_back(vector<int>(8,4));
    }
    stop_time= (int)clock();
    cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC) << endl;
    
    int x= 1;
}

/*-------- test custermized stxxl vector of vector --------*/
// an example from http://www.cplusplus.com/forum/general/268/
class myArray {
private:
    int size;
    int a[10];
public:
    int& operator[] (int x) {
        return a[x];
    }
    void print_array();   // I included this just to show the operator[] works!
};

void myArray::print_array()
{
    for (int j=0; j < 10; j++)
        cout << "array[" << j << "] = " << a[j] << "\n";
}

TEST_CASE("test custermized vector example from web", "[hide]")
{
    // create an instance of the myArray class
    myArray instance;
    
    // load instance.a[] with integers
    // NOTE: here we use instance[i] NOT instance.a[i]
    // instance.a[] wouldn't work as "int a[]" in the myArray class
    // is defined as PRIVATE!!!!!
    for (int i=0; i < 10; i++)
        instance[i] = i;
    
    // show that our operator worked by printing out the array values
    // using the myArray member function myArray::print_array()
    instance.print_array();
    
    cout << "\n\n";
    system("PAUSE");
}

// custermized stxxl vector of vector
// stxxl_v
/*class stxxl_v
{
public:
    stxxl_v(){}
    virtual ~stxxl_v(){}
    
public:
    size_t & operator [] (size_t j);
};

// stxxl_vv
class stxxl_vv
{
public:
    stxxl_vv(){}
    virtual ~stxxl_vv(){}
    
public:
    stxxl_v & operator [] (size_t i) { return this->dat_vec[i]; }
    void pileup_encode(string encode_file);
    void pileup_m5(string m5_file);
    
private:
    vector<stxxl_v> dat_vec;
    stxxl::vector<int64_t> dat;
    vector<int> pu_num;
    vector<int> pu_num_cum;
};

void stxxl_vv::pileup_encode(string encode_file)
{
    // get size pileup
    ifstream p_encode_file; open_infile(p_encode_file, encode_file);
    size_t pu_size = 0;
    while (true) {
        string buf;
        getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (auto i = 0; i < buf_vec.size(); ++i){
            if (buf_vec[i] > pu_size) pu_size = buf_vec[i];
        }
    }
    ++pu_size;
    p_encode_file.close();

    // get size of each locus
    vector<int> pu_num(pu_size, 0);
    size_t total_size = 0;
    open_infile(p_encode_file, encode_file);
    while (true) {
        string buf;
        getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (auto i = 0; i < buf_vec.size(); ++i){
            ++pu_num[buf_vec[i]];
        }
        total_size += buf_vec.size();
    }
    p_encode_file.close();
    this->pu_num = pu_num;
    
    this->dat.clear();
    this->dat.resize(total_size);
    
    vector<int> pu_num_cum(pu_num);
    for (auto i = 1; i < pu_num.size(); ++i)
        pu_num_cum[i] += pu_num_cum[i-1];
    this->pu_num_cum = pu_num_cum;
    int x = 1;
    
    // pileup on hard drive
    vector<int> pu_shift(pu_size, 0);
    open_infile(p_encode_file, encode_file);
    int64_t read_id = 0;
    while (true) {
        string buf;
        getline(p_encode_file, buf);
        if (p_encode_file.eof())
            break;
        vector<int> buf_vec = split_int(buf, '\t');
        for (auto i = 0; i < buf_vec.size(); ++i){
            int k = buf_vec[i] - 1;
            int cur_num_cum = k > 0 ? pu_num_cum[k] : 0;
            this->dat[cur_num_cum + pu_shift[buf_vec[i]]] = read_id;
            ++pu_shift[buf_vec[i]];
        }
        ++read_id;
    }
    p_encode_file.close();
    x = 1;
    
}

TEST_CASE("test custermized stxxl vector", "[hide]")
{
    string encode_file = "/Users/zhixingfeng/Dropbox/work/iGDA/paper/submission/nature_communication_revision/analysis/memory_reduce/pacbio_ecoli/results/igda/detect/qv0/realign.encode";
    string m5_file = "/Users/zhixingfeng/Dropbox/work/iGDA/paper/submission/nature_communication_revision/analysis/memory_reduce/pacbio_ecoli/results/igda/detect/qv0/realign.m5";
    
    stxxl_vv pu_var;
    stxxl_vv pu_read;
    
    pu_var.pileup_encode(encode_file);
    
}*/
