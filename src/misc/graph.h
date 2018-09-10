//
//  graph.h
//  iGDA
//
//  Created by Zhixing Feng on 9/10/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#ifndef iGDA_graph_h
#define iGDA_graph_h
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/transitive_reduction.hpp>
#include <boost/graph/graph_utility.hpp> // dumping graphs


using namespace boost;
struct IQsNode { };
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, IQsNode*> Graph;


inline void igda_transitive_reduction(const Graph in_g, Graph &out_g)
{
    std::map<Graph::vertex_descriptor, Graph::vertex_descriptor> g_to_tr;
    std::vector<size_t> id_map(num_vertices(in_g));
    std::iota(id_map.begin(), id_map.end(), 0u);
    
    transitive_reduction(in_g, out_g, make_assoc_property_map(g_to_tr), id_map.data());
}



#endif
