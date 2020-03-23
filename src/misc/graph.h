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

/*------------- BGL graph algorithm -------------*/
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
//void igda_transitive_reduction(const Graph in_g, Graph &out_g);

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

// get unambigious paths (no split)
set<GraphPath> get_unambigious_paths(const Graph &gp);

/*------------- igda graph algorithm --------------*/
// define igda graph type
struct IGDA_Vertex
{
    int64_t id; // ID the of vertex (0-based)
    int64_t start_locus; // start locus of the contig corresponding to the vertex
    int64_t end_locus; // end locus of the contig corresponding to the vertex
    int64_t npaths; // number of paths connected to the vertex
};
struct IGDA_Graph
{
    // adjacent matrix (out nodes), vertex ID starts from 0;
    map<int64_t, vector<IGDA_Vertex> > adj_mat;
    map<int64_t, vector<IGDA_Vertex> > adj_mat_in;
    
};

inline bool COMP_VERTEX_LESS (const IGDA_Vertex& a, const IGDA_Vertex& b)
{
    return (a.end_locus < b.end_locus);
}

// read igda graph from dot file and ann file
void load_igda_graph_from_file(IGDA_Graph &gp, string dot_file, string ann_file, bool is_sort = true);

// save igda graph to dot file
void save_igda_graph_to_file(const IGDA_Graph &gp, string dot_file);

// Eugene W. Myers's linear time complexity (Eugene W. Myers, The fragment assembly string graph, 2005)
void igda_tred(const IGDA_Graph &gp, IGDA_Graph &gp_tred);

// igda graph to boost graph
Graph convert_igda_graph_to_boost_graph(const IGDA_Graph &gp);

// get asscessible vertices from a vertex
void get_accessible_vertices(const IGDA_Graph &gp, unordered_set<int64_t> &accessible_vertices, int64_t start_vertex_id);

// get ambiguous paths (no more than two split)
set<vector<int64_t> > get_unambigious_paths_ms_core_legacy(const IGDA_Graph &gp, int64_t start_vertex_id, set<int64_t> &end_vertex_id);
set<vector<int64_t> > get_unambigious_paths_ms_legacy(const IGDA_Graph &gp);

void get_unambigious_paths_ms_core(const IGDA_Graph &gp, const unordered_set<int64_t> &accessible_vertices, int64_t vertex_id, vector<int64_t> &path, set<vector<int64_t> > &path_all,
                                   vector<bool> &visited, int &n_split, set<int64_t> &end_vertex_id);

#endif
