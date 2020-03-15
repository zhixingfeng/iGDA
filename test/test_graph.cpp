//
//  test_graph.cpp
//  iGDA
//
//  Created by Zhixing Feng on 10/29/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/misc.h"
#include "../src/misc/graph.h"
#include "../src/modules/assemble/assembler.h"

TEST_CASE("test accessing graph vertex set", "[hide]")
{
    /*--------------- Generate a graph ----------------*/
    //typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;
    
    // Make convenient labels for the vertices
    enum { A, B, C, D, E, N };
    const int num_vertices = N;
    const char* name = "ABCDE";
    
    // writing out the edges in the graph
    typedef std::pair<int, int> Edge;
    Edge edge_array[] =
    { Edge(A,B), Edge(A,D), Edge(C,A), Edge(D,C),
        Edge(C,E), Edge(B,D), Edge(D,E) };
    const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);
    
    // declare a graph object
    Graph g(num_vertices);
    
    // add the edges to the graph object
    for (int i = 0; i < num_edges; ++i)
        add_edge(edge_array[i].first, edge_array[i].second, g);
    
    cout << "graph = " << endl;
    boost::print_graph(g);
    cout << endl;
    
    /*---------------- Access vertex ------------------*/
    IndexMap index = get(boost::vertex_index, g);
    
    std::cout << "vertices(g) = ";

    pair<VertexIter, VertexIter> vp;
    for (vp = boost::vertices(g); vp.first != vp.second; ++vp.first){
        Vertex v = *vp.first;
        cout << index(v) << " ";
    }
    cout << endl << endl;
    
    /*----------------- Access edges ------------------*/
    EdgeIter ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei){
        cout << index(source(*ei, g)) << "->" << index(target(*ei, g)) << endl;
    }
    cout << endl;
    
    /*-------------- out edges ----------------*/
    
    cout << "out edges " << endl;
    for (vp = boost::vertices(g); vp.first != vp.second; ++vp.first){
        Vertex v = *vp.first;
        cout << index(v) << endl;
        OutEdgeIter out_i, out_end;
        for (boost::tie(out_i, out_end) = out_edges(v, g); out_i != out_end; ++out_i){
            Vertex cur_src = source(*out_i, g);
            Vertex cur_targ = target(*out_i, g);
            cout << index(cur_src) << "->" << index(cur_targ) << endl;
        }
    }
    
    /*------------- get vertices with no in-edge ---------------*/
    cout << "vertices with no out-edges are" << endl;
    for (vp = boost::vertices(g); vp.first != vp.second; ++vp.first){
        Vertex v = *vp.first;
        //cout << index(v) << endl;
        boost::graph_traits<Graph>::out_edge_iterator out_i, out_end;
        boost::tie(out_i, out_end) = out_edges(v, g);
        if (out_i == out_end)
            cout << index(v) << endl;
    }
    
}

TEST_CASE("test get_out_vertex", "[hide]")
{
    Graph gp;
    read_dot_file(gp, "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft.dot");
    pair<VertexIter, VertexIter> v_iter;
    IndexMap index = get(boost::vertex_index, gp);
    for (v_iter = boost::vertices(gp); v_iter.first != v_iter.second; ++v_iter.first){
        Vertex v = *v_iter.first;
        vector<Vertex> out_vertex = get_out_vertex(gp, v);
        cout << index(v) << " -> ";
        for (auto i = 0; i < out_vertex.size(); ++i){
            cout << index(out_vertex[i]) << ",";
        }
        cout << endl;
    }

    //print_graph(gp);
}



TEST_CASE("test get_vertices_no_inedge", "[hide]")
{
    Assembler assembler;
    Graph gp;
    assembler.ann_to_graph(gp, "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft");
    
    IndexMap index = get(boost::vertex_index, gp);
    vector<Vertex> v_no_inedge = get_vertices_no_inedge(gp);
    for (auto i = 0; i < v_no_inedge.size(); ++i)
        cout << index(v_no_inedge[i]) << endl;
    
}

TEST_CASE("test get_unambigious_paths", "[hide]")
{
    Graph gp;
    //read_dot_file(gp, "../results/pt_ann_assign_reads/align_to_consensus_trim.recode.ann.tested.ft.transitive_reduction.dot");
    read_dot_file(gp, "../results/sa/ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.encode.rdim.ann.tested.ft.tred.dot");
    
    
    //print_graph(gp);
    set<GraphPath> paths = get_unambigious_paths(gp);
    for (auto it = paths.begin(); it != paths.end(); ++it){
        for (auto i = 0; i < it->size(); ++i){
            cout << (*it)[i] << ',';
        }
        cout << endl;
    }
}

TEST_CASE("test load_igda_graph_from_file", "[hide]")
{
    string dot_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_500.dot";
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_500";
    
    string dot_file_rev = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_500.rev.dot";
    string ann_file_rev = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_500.rev";
    
    IGDA_Graph gp;
    load_igda_graph_from_file(gp, dot_file, ann_file);
    
    IGDA_Graph gp_rev;
    load_igda_graph_from_file(gp_rev, dot_file_rev, ann_file_rev);
    
    int x = 1;
    //read_dot_file(gp, dot_file);
}

TEST_CASE("test igda_tred (transitive reduction)", "[hide]")
{
    string dot_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_500.dot";
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_500";
    
    IGDA_Graph gp;
    load_igda_graph_from_file(gp, dot_file, ann_file);
    
    IGDA_Graph gp_tred;
    igda_tred(gp, gp_tred);
    
    save_igda_graph_to_file(gp_tred, dot_file + ".igda_tred.dot");
    
    Graph gp_bgl = convert_igda_graph_to_boost_graph(gp_tred);
    
    ofstream fs_graph(dot_file + ".igda_tred.blg.dot");
    boost::write_graphviz(fs_graph, gp_bgl);
    fs_graph.close();
    
    int x = 1;
}

TEST_CASE("test igda assemble positive dead loop issue", "[hide]")
{
    string dot_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_5000.tred.dot";
    string ann_file = "/Users/zhixingfeng/Dropbox/work/iGDA/development/test/test_tred/data/realign.ann.tested.ft.count.ft.head_5000";
    
    Assembler assembler;
    Graph gp;
    assembler.read_ann_results(ann_file);
    read_dot_file(gp, dot_file);
    assembler.assemble(gp, ann_file + ".assembled");
    
}







