#pragma once

#include <vector>
// reverse
#include <algorithm>

// to compute norms
#include "Eigen/Eigen"

#include "scomplex/Simplicial_Complex.h"

// graph libraries
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/property_map/vector_property_map.hpp"

typedef typename boost::adjacency_list<             //
    boost::vecS,                                    //
    boost::vecS,                                    //
    boost::undirectedS,                             //
    boost::no_property,                             // vertex property
    boost::property<boost::edge_weight_t, double>>  //
    Graph;

typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;

typedef std::vector<double> point_t;

Graph calculate_one_skelleton_graph(
    simplicial::SimplicialComplex<point_t>& scomp  //
    ) {
    int num_points{scomp.get_dimension(0)};
    int num_edges{scomp.get_dimension(1)};

    std::vector<std::pair<int, int>> edges(num_edges);
    std::vector<double> weights(num_edges);

    for (std::vector<int> edge : scomp.get_level(1)) {
        int src_index = edge.at(0);
        int trg_index = edge.at(1);

        point_t src_point = scomp.points.at(edge.at(0));
        point_t trg_point = scomp.points.at(edge.at(1));

        edges.push_back(std::make_pair(src_index, trg_index));

        Eigen::VectorXd diff(src_point.size());

        for (int i = 0; i < src_point.size(); i++)
            diff(i) = src_point.at(i) - trg_point.at(i);

        weights.push_back(diff.norm());
    }

    const Graph one_skelleton(edges.begin(),    //
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
std::vector<vertex_t> shortest_path(const Graph& g, vertex_t s, vertex_t t) {
    // make a predecessor map
    boost::vector_property_map<vertex_t> predecessors(boost::num_vertices(g));
    boost::dijkstra_shortest_paths(g, s,
                                   boost::predecessor_map(&predecessors[0]));
    // construct the path (going backwards from the target to the source)
    vertex_t it = t;
    std::vector<vertex_t> s_t_path{};
    // WARNING: the push_back and reverse strategy may be more efficient
    //
    // use a do-while block, we want the guard to be checked afterwards
    do {
        // push to the front
        s_t_path.insert(s_t_path.begin(), it);
        it = predecessors[it];
    } while (it != s);
    return s_t_path;
}

std::vector<vertex_t> complete_path(const Graph& g, std::vector<vertex_t> vec) {
    std::vector<vertex_t>::iterator s{vec.begin()};
    std::vector<vertex_t>::iterator t = std::next(s);

    // start calculatinng the full path
    std::vector<vertex_t> full_path;
    full_path.push_back(*s);
    while (t != vec.end()) {
        // get the next vertex
        std::cout << "looking at edge" << *s << ", " << *t << "\n";

        // calculate and add the path between s and t (contains vertices other
        // than s) to the full path
        std::vector<vertex_t> path_portion{shortest_path(g, *s, *t)};
        for (vertex_t vert : path_portion) {
            full_path.push_back(vert);
        }
        std::cout << "\n";
        std::advance(s, 1);
        std::advance(t, 1);
    }
    return full_path;
}
