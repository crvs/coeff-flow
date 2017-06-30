#include <iostream>

#include "scomplex/types.hpp"
#include "scomplex/graph_utils.hpp"

#include <vector>

using namespace std;


int main() {
    vector<pair<size_t,size_t>> graph_edges = {{0,1},{1,2},{2,3},{3,4},{0,4}};
    vector<double> graph_weights = {0,0,0,0,1};
    gsimp::graph_t g(graph_edges.begin(),graph_edges.end(),graph_weights.begin(),3);

    vector<size_t> p;
    p = gsimp::shortest_path(g,0,4);

    for(auto v : p ) cout << v << ' ';
    cout << '\n';

    return 0;
}
