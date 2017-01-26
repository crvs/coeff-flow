#ifndef UTILS
#define UTILS

#include <Eigen/Eigen>

#include <scomplex/Simplicial_Complex.h>

// graph libraries
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"

using namespace boost;

typedef typename property<edge_weight_t, double> Weight;
typedef
    typename adjacency_list<vecS, vecS, undirectedS, no_property, Weight> Graph;
typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

typedef std::vector<double> point_t;
Graph calculate_one_skelleton_graph(
    simplicial::SimplicialComplex<point_t>& scomp) {
    Graph one_skelleton(scomp.get_dimension(0));

    for (std::vector<int> edge : scomp.get_level(1)) {
        point_t src = scomp.points.at(edge.at(0));
        point_t trg = scomp.points.at(edge.at(1));

        Eigen::VectorXd diff(src.size());
        for (int i = 0; i < src.size(); i++) {
            diff(i) = src.at(i) - trg.at(i);
        }

        double norm = diff.norm();
        add_edge(edge.at(0), edge.at(1), Weight(norm), one_skelleton);
    }
    return one_skelleton;
}

// calculate shortest paths:
auto shortest_path(Graph g, Vertex s, Vertex t) {
    std::cout << "I'm calculating shortest paths\n";
    dijkstra_shortest_paths(Graph g, Vertex v);
    return 0;
}

#endif
