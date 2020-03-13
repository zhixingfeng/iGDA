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
#include "./basic.h"

//using namespace boost;
struct IQsNode { };
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, IQsNode*> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIter;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIter;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
typedef boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;
typedef vector<Vertex> GraphPath;

// transtive reduction of BGL is INCORRECT, NEVER USE IT!!!
void igda_transitive_reduction(const Graph in_g, Graph &out_g);

// read dot file
void read_dot_file(Graph &gp, string dot_file);

// get number of vertices
size_t get_num_vertices(const Graph &gp);

// get vertices with no in-edge
vector<Vertex> get_vertices_no_inedge(const Graph &gp);

// get number of out-edges for a vertex
vector<Vertex> get_out_vertex(const Graph &gp, const Vertex &v);

// travel through a path until hit a vertex with no or more than 1 out edges
GraphPath travel_path(const Graph &gp, const Vertex &v_start, int i = 0);

// get unambigious paths
set<GraphPath> get_unambigious_paths(const Graph &gp);

/*------------- define igda graph --------------*/
// define igda graph type
struct IGDA_Vertex
{
    int64_t id;
    int64_t start_locus;
    int64_t end_locus;
};
struct IGDA_Graph
{
    // adjacent matrix, vertex ID starts from 0;
    map<int64_t, vector<IGDA_Vertex> > adj_mat;
};

/*struct COMP_VERTEX_LESS
{
    inline bool operator() (const IGDA_Vertex& a, const IGDA_Vertex& b)
    {
        return (a.end_locus < b.end_locus);
    }
};*/
inline bool COMP_VERTEX_LESS (const IGDA_Vertex& a, const IGDA_Vertex& b)
{
    return (a.end_locus < b.end_locus);
}

// read igda graph from dot file and ann file
void load_igda_graph_from_file(IGDA_Graph &gp, string dot_file, string ann_file, bool is_sort = true);



#endif
