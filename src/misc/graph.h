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
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>

#include "../../include/headers.h"

//using namespace boost;
struct IQsNode { };
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, IQsNode*> Graph;

void igda_transitive_reduction(const Graph in_g, Graph &out_g);


#endif
