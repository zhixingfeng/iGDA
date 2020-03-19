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
/*void igda_transitive_reduction(const Graph in_g, Graph &out_g)
{
    std::map<Graph::vertex_descriptor, Graph::vertex_descriptor> g_to_tr;
    std::vector<size_t> id_map(num_vertices(in_g));
    std::iota(id_map.begin(), id_map.end(), 0u);
    
    transitive_reduction(in_g, out_g, boost::make_assoc_property_map(g_to_tr), id_map.data());
}*/

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
    unordered_set<Vertex> checked_vertice;
    while(v_active.size()>0){
        
        Vertex cur_vertex = v_active.top();
        
        // debug start
        //cout << "size of v_active = " << v_active.size() << endl;
        //cout << "current v_active = " << v_active.top() << endl << endl;
        //getchar();
        // debug end
        
        auto it_find_vertex = checked_vertice.find(cur_vertex);
        if (it_find_vertex == checked_vertice.end()){
            checked_vertice.insert(cur_vertex);
        }else{
            //cout << "cur_vertex " << cur_vertex << " exists" << endl;
            //getchar();
            v_active.pop();
            continue;
        }

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
            // debug start
            //for (auto t : cur_path)
            //    cout << t << " ";
            //cout << endl << endl;
            //getchar();
            // debug end
            paths.insert(cur_path);
            size_t num_out_vertex = get_out_vertex(gp, cur_path.back()).size();
            
            if (num_out_vertex == 1)
                throw runtime_error("get_unambigious_paths(): num_out_vertex == 1");
            
            if (cur_path.back() != cur_vertex && get_out_vertex(gp, cur_path.back()).size() > 1 &&
                checked_vertice.find(cur_path.back()) == checked_vertice.end())
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
        gp.adj_mat_in[i] = vector<IGDA_Vertex>();
    }
    for (auto i = 0; i < edges.size(); ++i){
        // record adj_mat
        IGDA_Vertex cur_vertex;
        cur_vertex.id = edges[i].second;
        cur_vertex.start_locus = ann_data[edges[i].second].start;
        cur_vertex.end_locus = ann_data[edges[i].second].end;
        gp.adj_mat[edges[i].first].push_back(cur_vertex);
        
        // record adj_mat_in
        IGDA_Vertex cur_vertex_in;
        cur_vertex_in.id = edges[i].first;
        cur_vertex_in.start_locus = ann_data[edges[i].first].start;
        cur_vertex_in.end_locus = ann_data[edges[i].first].end;
        gp.adj_mat_in[edges[i].second].push_back(cur_vertex_in);
    }
    
    // sort linked vertex
    if (is_sort){
        for (auto it = gp.adj_mat.begin(); it != gp.adj_mat.end(); it++){
            std::stable_sort(it->second.begin(), it->second.end(), COMP_VERTEX_LESS);
        }
    }
}
// save igda graph to dot file
void save_igda_graph_to_file(const IGDA_Graph &gp, string dot_file)
{
    ofstream fs_dot_file; open_outfile(fs_dot_file, dot_file);
    fs_dot_file << "digraph G {" << endl;
    for (auto it = gp.adj_mat.begin(); it != gp.adj_mat.end(); ++it){
        
        if (it->second.size() == 0){
            fs_dot_file <<  '\t' << it->first << ";" << endl;
        }else{
            for (auto i = 0; i < it->second.size(); ++i){
                fs_dot_file << '\t' << it->first << " -> " << it->second[i].id << ";" << endl;
            }
        }
    }
    fs_dot_file << "}";
    fs_dot_file.close();
    
}

// Eugene W. Myers's linear time complexity transitive reduction algorithm (Eugene W. Myers, The fragment assembly string graph, 2005)
void igda_tred(const IGDA_Graph &gp, IGDA_Graph &gp_tred)
{
    //gp_tred = gp;
    // initialize reduced graph
    for (auto it = gp.adj_mat.begin(); it != gp.adj_mat.end(); ++it)
        gp_tred.adj_mat[it->first] = vector<IGDA_Vertex>();
    
    // initialize the status of each vertex
    enum VertexStatus {vacant, inplay, eliminated};
    vector<VertexStatus> status(gp.adj_mat.size(), vacant);
    
    // for each vertex 
    for (auto it_i = gp.adj_mat.begin(); it_i != gp.adj_mat.end(); ++it_i){
        // if out degree of the current vertex is 0, skip
        if (it_i->second.size() == 0) continue;
        
        // mark the directed vertex of the current vertex as inplay
        for (auto i = 0; i < it_i->second.size(); ++i)
            status[it_i->second[i].id] = inplay;
        
        // get the right most locus
        int64_t longest = it_i->second[it_i->second.size() - 1].end_locus;
        
        // find reducible directed vertex
        //cout << it_i->first << endl;
        for (auto i = 0; i < it_i->second.size(); ++i){
            int64_t cur_vertex = it_i->second[i].id;
            auto it_j = gp.adj_mat.find(cur_vertex);
            if (it_j == gp.adj_mat.end())
                throw runtime_error("igda_tred(): it_j == gp.adj_mat.end()");
            for (auto j = 0; j < it_j->second.size(); ++j){
                if (it_j->second[j].end_locus > longest) break;
                if (it_j->second[j].id >= status.size()){
                    cout << "fault:" << endl;
                    cout << it_j->first << endl;
                    cout << it_j->second[j].id << endl;
                    throw runtime_error("igda_tred(): it_j->second[j].id >= status.size()");
                }
                if (status[it_j->second[j].id] == inplay){
                    status[it_j->second[j].id] = eliminated;
                }
            }
        }
        
        // remove reducible edges
        for (auto i = 0; i < it_i->second.size(); ++i){
            if (status[it_i->second[i].id] != eliminated){
                gp_tred.adj_mat[it_i->first].push_back(it_i->second[i]);
            }
            status[it_i->second[i].id] = vacant;
        }
    }
    
}

Graph convert_igda_graph_to_boost_graph(const IGDA_Graph &gp)
{
    Graph gp_bgl;
    
    for (auto it = gp.adj_mat.begin(); it != gp.adj_mat.end(); ++it){
        boost::add_vertex(gp_bgl);
    }
  
    for (auto it = gp.adj_mat.begin(); it != gp.adj_mat.end(); ++it){
        for (auto i = 0; i < it->second.size(); ++i){
            boost::add_edge(it->first, it->second[i].id, gp_bgl);
        }
    }
    return gp_bgl;
}

void get_accessible_vertices(const IGDA_Graph &gp, unordered_set<int64_t> &accessible_vertices, int64_t start_vertex_id)
{
    // init
    stack<int64_t> v_active;
    v_active.push(start_vertex_id);
    vector<int64_t> v_visited(gp.adj_mat.size(), false);
    
    // depth first search
    while(!v_active.empty()){
        cout << v_active.size() << ":" << v_active.top() << endl;
        
        int64_t cur_v = v_active.top();
        v_active.pop();
        
        if (!v_visited[cur_v]){
            accessible_vertices.insert(cur_v);
            v_visited[cur_v] = true;
        }
        
        auto it = gp.adj_mat.find(cur_v);
        if (it == gp.adj_mat.end())
            throw runtime_error("get_npaths_between_vertices_core(): fail to find cur_v in graph");
        
        vector<IGDA_Vertex> cur_adj_mat = it->second;
        for (auto i = 0; i < cur_adj_mat.size(); ++i){
            if (!v_visited[cur_adj_mat[i].id])
                v_active.push(cur_adj_mat[i].id);
        }
    }
}
