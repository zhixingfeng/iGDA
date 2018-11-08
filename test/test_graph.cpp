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




TEST_CASE("test transitive reduction", "[hide]")
{
    Graph g;
    enum {a,b,c,d,e};
    add_edge(a,b,g);
    add_edge(a,c,g);
    add_edge(a,d,g);
    add_edge(a,e,g);
    add_edge(b,d,g);
    add_edge(c,d,g);
    add_edge(c,e,g);
    add_edge(d,e,g);
    
    Graph tr;
    igda_transitive_reduction(g, tr);
    /*std::map<Graph::vertex_descriptor, Graph::vertex_descriptor> g_to_tr;
     std::vector<size_t> id_map(num_vertices(g));
     std::iota(id_map.begin(), id_map.end(), 0u);
     
     transitive_reduction(g, tr, make_assoc_property_map(g_to_tr), id_map.data());*/
    
    boost::print_graph(g);
    std::cout << "----------------------------\n";
    boost::print_graph(tr);
    
    /*// generating graphviz files
     { std::ofstream dot("g.dot");  write_graphviz(dot, g); }
     { std::ofstream dot("tr.dot"); write_graphviz(dot, tr); }*/
}

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








