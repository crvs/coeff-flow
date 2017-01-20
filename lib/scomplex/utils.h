#ifndef UTILS
#define UTILS

#include <Eigen/Eigen>
#include <scomplex/Simplicial_Complex.h>

// graph libraries
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"

typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
    boost::property<boost::edge_weight_t, double>> Graph;

template <typename point_t>
Graph get_one_skelleton_graph(simplicial::SimplicialComplex<point_t>& scomp) {
    Graph one_skelleton(scomp.get_dimension(0));

    for (std::vector<int> edge : scomp.get_level(1)) {
        std::cout << "adding edges: " << std::endl;
        for (auto c : edge) {
            std::cout << c << " ";
        }

        for (auto c : scomp.points.at(edge.at(0))) {
            std::cout << c << " ";
        }
        std::cout << std::endl;

        point_t src = scomp.points.at(edge.at(0));
        /*
        Eigen::VectorXd src_v(scomp.points.at(edge.at(0)).data());

        Eigen::VectorXd trg(scomp.points.at(edge.at(1)).data());
        Eigen::VectorXd diff = src - trg;
        double norm = .1;  // diff.norm();

        auto diff = src - trg;
        double norm = diff.dot(diff);
        norm = std::sqrt(norm);
        boost::add_edge(edge.at(0), edge.at(1), norm, one_skelleton);
        */
    }
    return one_skelleton;
}

#endif
