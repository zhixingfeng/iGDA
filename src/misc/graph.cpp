//
//  graph.cpp
//  iGDA
//
//  Created by Zhixing Feng on 10/29/18.
//  Copyright (c) 2018 Zhixing Feng. All rights reserved.
//

#include "./graph.h"
#include "../modules/assemble/assembler.h"

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

size_t get_num_vertices(const Graph &gp)
{
    size_t num_vertices = 0;
    pair<VertexIter, VertexIter> v_iter;
    for (v_iter = boost::vertices(gp); v_iter.first != v_iter.second; ++v_iter.first)
        ++num_vertices;
    return num_vertices;
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


vector<Vertex> get_out_vertex(const Graph &gp, const Vertex &v)
{
    OutEdgeIter out_edge_i, out_edge_end;
    vector<Vertex> out_vertext;
    for (boost::tie(out_edge_i, out_edge_end) = out_edges(v, gp); out_edge_i != out_edge_end; ++out_edge_i)
        out_vertext.push_back(target(*out_edge_i, gp));
    return out_vertext;
}


GraphPath travel_path(const Graph &gp, const Vertex &v_start, int i)
{
    GraphPath path;
    Vertex cur_v = v_start;
    while(true){
        path.push_back(cur_v);
        vector<Vertex> out_vertex = get_out_vertex(gp, cur_v);
        
        if (out_vertex.size() == 0)
            break;
        
        if (out_vertex.size() > 1){
            if (cur_v != v_start)
                break;
            else{
                if (i >= out_vertex.size())
                    throw runtime_error("travel_path(): i >= out_vertex.size()");
                cur_v = out_vertex[i];
            }
        }else{
            cur_v = out_vertex[0];
        }
        
        
    }
    return path;
}


set<GraphPath> get_unambigious_paths(const Graph &gp)
{
    //IndexMap index = get(boost::vertex_index, gp);
    
    set<GraphPath> paths;
    vector<Vertex> v_no_inedge = get_vertices_no_inedge(gp);
    stack<Vertex, vector<Vertex> > v_active (v_no_inedge);
    
    // if find the second break points in the paths, stop, record the paths and add the break point to v_active
    while(v_active.size()>0){
        Vertex cur_vertex = v_active.top();
        v_active.pop();
        
        vector<Vertex> cur_out_vertex = get_out_vertex(gp, cur_vertex);
        
        if (cur_out_vertex.size() == 0){
            GraphPath cur_path;
            cur_path.push_back(cur_vertex);
            paths.insert(cur_path);
            continue;
        }
        
        for (auto i = 0; i < cur_out_vertex.size(); ++i){
            GraphPath cur_path = travel_path(gp, cur_vertex, i);
            paths.insert(cur_path);
            size_t num_out_vertex = get_out_vertex(gp, cur_path.back()).size();
            
            if (num_out_vertex == 1)
                throw runtime_error("get_unambigious_paths(): num_out_vertex == 1");
            
            if (cur_path.back() != cur_vertex && get_out_vertex(gp, cur_path.back()).size() > 1)
                v_active.push(cur_path.back());
        }
        
    }
    
    return paths;
}

// read igda graph from dot file and ann file
void load_igda_graph_from_file(IGDA_Graph &gp, string dot_file, string ann_file, bool is_sort)
{
    // load ann file
    Assembler assembler;
    assembler.read_ann_results(ann_file);
    vector<ConsensusSeq> ann_data = assembler.get_ann_result();
    
    // load graph from dot file
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
    
    // construct igda graph
    for (auto i = 0; i <= node_max; ++i){
        gp.adj_mat[i] = vector<IGDA_Vertex>();
    }
    for (auto i = 0; i < edges.size(); ++i){
        IGDA_Vertex cur_vertex;
        cur_vertex.id = edges[i].second;
        cur_vertex.start_locus = ann_data[edges[i].second].start;
        cur_vertex.end_locus = ann_data[edges[i].second].end;
        gp.adj_mat[edges[i].first].push_back(cur_vertex);
    }
    
    // sort linked vertex
    if (is_sort){
        for (auto it = gp.adj_mat.begin(); it != gp.adj_mat.end(); it++){
            std::stable_sort(it->second.begin(), it->second.end(), COMP_VERTEX_LESS);
        }
    }
}

