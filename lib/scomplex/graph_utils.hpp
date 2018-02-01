#pragma once

#include <iostream>
#include <vector>
#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>

// graph libraries
#include "boost/config.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/property_map/vector_property_map.hpp"

namespace gsimp {

using namespace boost;
/*
typedefs:
    graph_t
functions:
    calculate_one_skelleton_graph(simplicial_complex& s_comp )
    complete_path(const graph_t& g, std::vector<vertex_t> vec)
    shortest_path(const graph_t& g, vertex_t s, vertex_t t)
*/

typedef adjacency_list<               //
    vecS,                             //
    vecS,                             //
    undirectedS,                      //
    no_property,                      // vertex property
    property<edge_weight_t, double>>  //
    graph_t;

typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
typedef std::pair<size_t, size_t> Edge;

graph_t calculate_one_skelleton_graph(simplicial_complex& s_comp  //
                                      ) {
    size_t num_points(s_comp.get_level_size(0));
    size_t num_edges(s_comp.get_level_size(1));

    std::vector<Edge> g_edges(num_edges);
    std::vector<double> weights(num_edges);

    auto edges = s_comp.get_level(1);

    graph_t one_skelleton(num_points);

    for (size_t ind = 0; ind < edges.size(); ++ind) {
        size_t src_index = edges.at(ind).at(0);
        size_t trg_index = edges.at(ind).at(1);

        point_t src_point = s_comp.get_point(src_index);
        point_t trg_point = s_comp.get_point(trg_index);

        g_edges.push_back(std::make_pair(src_index, trg_index));

        double norm = 0;
        for (int i = 0; i < src_point.size(); ++i) { 
            double sq_d = src_point.at(i) - trg_point.at(i);
            sq_d = pow(sq_d,2);
            norm += sq_d;
        }
        norm = sqrt(norm); //
        // (not really) TODO (but we want to check this)
        // penalize the use of lots of edges
        boost::add_edge(src_index,trg_index,norm,one_skelleton);
    }


   /* graph_t one_skelleton(g_edges.begin(),  //
    *                             g_edges.end(),    //
    *                             weights.begin(),  //
    *                             num_points);      //
    *
    */
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
std::vector<vertex_t> shortest_path(graph_t& g, vertex_t s, vertex_t t) {
    // make a predecessor map
    // property_map<graph_t, edge_weight_t>::type weight_map = get(edge_weight, g);
    std::vector<vertex_t> predecessors(num_vertices(g));
    std::vector<double> distances(num_vertices(g));
    dijkstra_shortest_paths(g, s, predecessor_map(&predecessors[0]).distance_map(&distances[0]));
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

std::vector<vertex_t> complete_path(graph_t& g,
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
        full_path.insert(full_path.end(), path_portion.begin(),
                         path_portion.end());
        std::advance(s, 1);
        std::advance(t, 1);
    }
    return full_path;
}
};
