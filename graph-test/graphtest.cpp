#include <iostream>

#include "boost/array.hpp"
#include "boost/graph/visitors.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/named_function_params.hpp"
#include "boost/graph/breadth_first_search.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"

int main() {
    typedef typename boost::adjacency_list<  //
        boost::setS,                         //
        boost::vecS,                         //
        boost::undirectedS,                  //
        boost::no_property,                  //
        boost::property<boost::edge_weight_t, double>> graph_t;

    graph_t g;

    for (int v = 0; v < 4; v++) {
        boost::add_vertex(g);
    }

    /*
    std::pair<                                     //
        boost::adjacency_list<>::vertex_iterator,  //
        boost::adjacency_list<>::vertex_iterator> verts = boost::vertices(g);
    */
    boost::add_edge(0, 1, 1.3, g);
    boost::add_edge(1, 3, 1.5, g);
    boost::add_edge(2, 3, 1.6, g);
    std::cout << "\nnumber of vertices: "  //
              << boost::num_vertices(g)    //
              << "\nnumber of edges: "     //
              << boost::num_edges(g)       //
              << "\nvertices: \n";

    typename boost::adjacency_list<>::vertex_iterator it, end;
    std::tie(it, end) = boost::vertices(g);
    std::copy(it, end,
              std::ostream_iterator<boost::adjacency_list<>::vertex_descriptor>{
                  std::cout, ", "});
    std::cout << "\n edges: \n";
    graph_t::edge_iterator ies, ees;
    std::tie(ies, ees) = boost::edges(g);
    std::copy(ies, ees,
              std::ostream_iterator<graph_t::edge_descriptor>{std::cout, "\n"});

    boost::array<int, 4> distances{{0}};

    boost::array<int, 4> predecessors;
    predecessors[1] = 1;
    boost::breadth_first_search(                 //
        g,                                       //
        1,                                       //
        boost::visitor(                          //
            boost::make_bfs_visitor(             // visitor:
                std::make_pair(                  // doing both at the same time
                    boost::record_distances(     // distance recorder
                        distances.begin(),       // distance array
                        boost::on_tree_edge{}),  // event
                    boost::record_predecessors(  // predecessor recorder
                        predecessors.begin(),    // predecessor array
                        boost::on_tree_edge{}))  // event
                )));                             //

    std::copy(distances.begin(), distances.end(),
              std::ostream_iterator<double>{std::cout, "\n"});

    int p = 2;
    std::cout << "path from 2: ";
    while (p != 1) {
        std::cout << p << ", ";
        p = predecessors[p];
    }
    std::cout << p << '\n';

    boost::array<int, 4> directions;
    boost::dijkstra_shortest_paths(  //
        g,                           //
        1,                           //
        boost::predecessor_map(directions.begin()));

    p = 2;
    std::cout << "dijkstra path from 2: ";
    while (p != 1) {
        std::cout << p << ", ";
        p = directions[p];
    }
    std::cout << p << "\n";

    return 0;
}
