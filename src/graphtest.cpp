
#include <iostream>

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    Graph;

int main() {
    Graph g(10);
    boost::add_edge(0, 1, g);

    return 0;
}
