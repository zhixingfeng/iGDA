//
//  test_BGL.cpp
//  iGDA
//
//  Created by Zhixing Feng on 9/7/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "../include/catch.hpp"
#include "../include/headers.h"
#include "../src/misc/misc.h"
#include "../src/misc/graph.h"


TEST_CASE("test transitive reduction")
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
    
    print_graph(g);
    std::cout << "----------------------------\n";
    print_graph(tr);
    
    /*// generating graphviz files
    { std::ofstream dot("g.dot");  write_graphviz(dot, g); }
    { std::ofstream dot("tr.dot"); write_graphviz(dot, tr); }*/
}