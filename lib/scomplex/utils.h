#ifndef UTILS
#define UTILS

#include <vector>
#include <algorithm>

#include "Eigen/Eigen"

#include "scomplex/Simplicial_Complex.h"

// graph libraries
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/property_map/vector_property_map.hpp"

using namespace boost;

typedef adjacency_list<               //
    vecS,                             //
    vecS,                             //
    undirectedS,                      //
    no_property,                      //
    property<edge_weight_t, double>>  //
    Graph;

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

typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

// calculate shortest paths:
std::vector<Vertex> shortest_path(const Graph& g, Vertex s, Vertex t) {
    boost::vector_property_map<Vertex> predecessors(num_vertices(g));

    dijkstra_shortest_paths(g, s, predecessor_map(&predecessors[0]));

    Vertex it = t;
    std::vector<Vertex> s_t_path{};
    while (it != s) {
        s_t_path.push_back(it);
        it = predecessors[it];
    }
    s_t_path.push_back(it);
    std::reverse(s_t_path.begin(), s_t_path.end());
    return s_t_path;

    // std::cout << "I'm calculating shortest paths\n";
}

std::vector<Vertex> complete_path(const Graph& g, std::vector<Vertex> vec) {
    std::vector<Vertex>::iterator s{vec.begin()};
    std::list<Vertex> full_path{*s};
    while (s != vec.end()) {
        std::vector<Vertex>::iterator t = std::next(s);
        // calculate the path between s and t
        std::vector<Vertex> path_portion{shortest_path(g, *s, *t)};

        // add the path between s and t to the full path
        for (Vertex vert : path_portion) {
            // allways skip the first element since we are calculating the paths
            // between s,t then s',t' where s' = t and that would make use of a
            // non-existing edge.
            if (vert != *path_portion.begin()) {
                full_path.push_back(vert);
            }
        }

        std::advance(s, 1);
    }
    // initialize with size so we don't have to relocate
    std::vector<Vertex> full_path_vect(full_path.size());
    for (Vertex vert : full_path) {
        full_path_vect.push_back(vert);
    }

    return full_path_vect;
}

#endif
