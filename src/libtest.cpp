#include <scomplex/Simplicial_Complex.h>
#include <scomplex/utils.h>
#include <scomplex/trace.h>
#include <scomplex/qhull_parsing.hpp>
#include <scomplex/path_snapper.hpp>

#include "boost/graph/graph_traits.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <tuple>

typedef std::vector<double> point_t;
typedef std::list<int> simp_t;
typedef std::pair<int, std::vector<int>> chain_t;

int f(point_t a) { return a.at(0) < .5 ? 0 : 1; }

int main() {
    std::list<point_t> point_list;
    std::list<simp_t> simp_list;
    // ---

    std::vector<point_t> points_v;
    std::vector<cell_t> cells_v;

    std::string filename;
    std::cout << "insert file name: ";
    std::cin >> filename;
    std::tie(points_v, cells_v) =
        // parse_qhull_file("/home/crvs/dev-cpp/examples/qhull-test/qh-test.dat");
        parse_qhull_file(filename);

    simplicial::SimplicialComplex rsc(points_v, cells_v);
    snap::path_snapper snapper(rsc);
    // snap::path_snapper snapper(points_v, cells_v);

    std::cout << "\n";
    auto path = snapper.snap_path({{-1, -1}, {1, 0}, {1, 1}});
    std::copy(path.begin(), path.end(),
              std::ostream_iterator<size_t>(std::cout, " "));
    return 0;
}
