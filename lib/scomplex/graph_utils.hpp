#pragma once

#include <vector>
#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>

// graph libraries
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/property_map/vector_property_map.hpp"

namespace gsimp {

/*
typedefs:
    graph_t
functions:
    calculate_one_skelleton_graph(simplicial_complex& s_comp )
    complete_path(const graph_t& g, std::vector<vertex_t> vec)
    shortest_path(const graph_t& g, vertex_t s, vertex_t t)
*/

typedef typename boost::adjacency_list<             //
    boost::vecS,                                    //
    boost::vecS,                                    //
    boost::undirectedS,                             //
    boost::no_property,                             // vertex property
    boost::property<boost::edge_weight_t, double>>  //
    graph_t;

typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;

graph_t calculate_one_skelleton_graph(simplicial_complex& s_comp  //
                                      ) {
    size_t num_points(s_comp.get_level_size(0));
    size_t num_edges(s_comp.get_level_size(1));

    std::vector<std::pair<int, int>> edges(num_edges);
    std::vector<double> weights(num_edges);

    for (cell_t edge : s_comp.get_level(1)) {
        int src_index = edge.at(0);
        int trg_index = edge.at(1);

        point_t src_point = s_comp.get_points().at(edge.at(0));
        point_t trg_point = s_comp.get_points().at(edge.at(1));

        edges.push_back(std::make_pair(src_index, trg_index));

        double norm = 0;
        for (int i = 0; i < src_point.size(); ++i)
            norm += pow(src_point.at(i) - trg_point.at(i), 2);
        weights.push_back(sqrt(norm));
    }

    const graph_t one_skelleton(edges.begin(),    //
                                edges.end(),      //
                                weights.begin(),  //
                                num_points);      //
    return one_skelleton;
}

/**
* @brief  find the shortest path between two vertices in a graph.
*
* @param g: graph (passed by ref)
* @param s: source vertex (type: decltype(g)::vertex_descriptor)
* @param t: source vertex (type: decltype(g)::vertex_descriptor)
*
* @returns: a vector of vertices
*            (type: std::vector<decltype(g)::vertex_descriptor>)
*/
std::vector<vertex_t> shortest_path(const graph_t& g, vertex_t s, vertex_t t) {
    // make a predecessor map
    boost::vector_property_map<vertex_t> predecessors(boost::num_vertices(g));
    boost::dijkstra_shortest_paths(g, s,
                                   boost::predecessor_map(&predecessors[0]));
    // construct the path (going backwards from the target to the source)
    vertex_t it = t;
    std::vector<vertex_t> s_t_path{};
    do {
        s_t_path.push_back(it);
        it = predecessors[it];
    } while (it != s);
    std::reverse(s_t_path.begin(), s_t_path.end());
    return s_t_path;
}

std::vector<vertex_t> complete_path(const graph_t& g,
                                    std::vector<vertex_t> vec) {
    std::vector<vertex_t>::iterator s{vec.begin()};
    std::vector<vertex_t>::iterator t = std::next(s);

    // start calculatinng the full path
    std::vector<vertex_t> full_path;
    full_path.push_back(*s);
    while (t != vec.end()) {
        // get the next vertex calculate and add the path between s and t
        // (contains vertices other than s) to the full path
        std::vector<vertex_t> path_portion{shortest_path(g, *s, *t)};
        // quick vector merging
        full_path.insert(full_path.end(),path_portion.begin(),path_portion.end());
        std::advance(s, 1);
        std::advance(t, 1);
    }
    return full_path;
}
};
