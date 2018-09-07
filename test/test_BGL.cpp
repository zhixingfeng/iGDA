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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/transitive_reduction.hpp>
#include <iostream>

//#include <boost/graph/graph_utility.hpp> // dumping graphs
//#include <boost/graph/graphviz.hpp>      // generating pictures

using namespace boost;

struct IQsNode { };
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, IQsNode*> Graph;

Graph make_random();


TEST_CASE("test transitive reduction")
{
    //Graph const g = make_random();
    Graph const g;
    Graph tr;
    std::map<Graph::vertex_descriptor, Graph::vertex_descriptor> g_to_tr;
    std::vector<size_t> id_map(num_vertices(g));
    std::iota(id_map.begin(), id_map.end(), 0u);
    
    transitive_reduction(g, tr, make_assoc_property_map(g_to_tr), id_map.data());
    
    /*print_graph(g);
    std::cout << "----------------------------\n";
    for (auto& e : g_to_tr)
        std::cout << "Mapped " << e.first << " to " << e.second << "\n";
    std::cout << "----------------------------\n";
    print_graph(tr);
    
    // generating graphviz files
    { std::ofstream dot("g.dot");  write_graphviz(dot, g); }
    { std::ofstream dot("tr.dot"); write_graphviz(dot, tr); }*/
}