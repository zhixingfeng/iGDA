//
//  graph.cpp
//  iGDA
//
//  Created by Zhixing Feng on 10/29/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "./graph.h"

// transtive reduction of BGL is INCORRECT, NEVER USE IT!!!
void igda_transitive_reduction(const Graph in_g, Graph &out_g)
{
    std::map<Graph::vertex_descriptor, Graph::vertex_descriptor> g_to_tr;
    std::vector<size_t> id_map(num_vertices(in_g));
    std::iota(id_map.begin(), id_map.end(), 0u);
    
    transitive_reduction(in_g, out_g, boost::make_assoc_property_map(g_to_tr), id_map.data());
}

// load graph from dot file. IMPORTANT: gp should be empty!!
void read_dot_file(Graph &gp, string dot_file)
{
    typedef std::pair<int, int> Edge;
    vector<Edge> edges;
    int64_t node_max = -1;
    
    // read dot file
    ifstream fs_infile;
    open_infile(fs_infile, dot_file);
    
    int64_t n_line = 0;
    while(true){
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof())
            break;
        
        ++n_line;
        
        size_t pos_found = buf.find(";");
        if (pos_found == string::npos)
            continue;
        
        // read current node and edge
        buf.replace(pos_found, 1, "");
        
        stringstream ss_buf;
        ss_buf << buf;
        vector<string> buf_vec;
        
        while(!ss_buf.eof()){
            string cur_buf;
            ss_buf >> cur_buf;
            buf_vec.push_back(cur_buf);
        }
        if (buf_vec.size() != 1 && buf_vec.size() != 2 && buf_vec.size() != 3)
            throw runtime_error(":read_dot_file: buf_vec.size() != 1 && buf_vec.size() != 3");
        
        
        // update node if there is only one node in the current line
        int64_t out_node = stoll(buf_vec[0]);
        if (out_node > node_max)
            node_max = out_node;
        
        // add edge if there are two nodes in the current line
        int64_t in_node = -1;
        if (buf_vec.size() == 3){
            if (buf_vec[1] != "->")
                throw runtime_error("read_dot_file: buf_vec[1] != \"->\"");
            in_node = stoll(buf_vec[2]);
            
            if (in_node > node_max)
                node_max = in_node;
            
            edges.push_back(Edge(out_node, in_node));
        }
        
        // deal with output of write_graphviz of BGL
        if (buf_vec.size() == 2){
            size_t cur_pos = buf_vec[0].find("->");
            if (cur_pos != string::npos){
                buf_vec[0].replace(cur_pos, 2, "\t");
                vector<string> cur_buf_vec = split(buf_vec[0], '\t');
                if (cur_buf_vec.size() != 2)
                    throw runtime_error("read_dot_file: cur_buf_vec.size() != 2");
                
                out_node = stoll(cur_buf_vec[0]);
                in_node = stoll(cur_buf_vec[1]);
                
                if (out_node > node_max)
                    node_max = out_node;
                
                if (in_node > node_max)
                    node_max = in_node;
                
                edges.push_back(Edge(out_node, in_node));
            }
        }
    }
    
    fs_infile.close();
    
    // update graph
    for (auto i = 0; i <= node_max; ++i)
        boost::add_vertex(gp);
    
    for (auto i = 0; i < edges.size(); ++i){
        boost::add_edge(edges[i].first, edges[i].second, gp);
    }
    
}

vector<Vertex> get_vertices_no_inedge(const Graph &gp)
{
    vector<Vertex> v_no_inedge;
    pair<VertexIter, VertexIter> v_iter;
    for (v_iter = boost::vertices(gp); v_iter.first != v_iter.second; ++v_iter.first){
        Vertex v = *v_iter.first;
        InEdgeIter in_edge_i, in_edge_end;
        boost::tie(in_edge_i, in_edge_end) = in_edges(v, gp);
        if (in_edge_i == in_edge_end)
            v_no_inedge.push_back(v);
    }
    return v_no_inedge;
}





