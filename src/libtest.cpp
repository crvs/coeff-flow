#include <scomplex/Simplicial_Complex.h>
#include <scomplex/utils.h>
#include <scomplex/trace.h>
#include <scomplex/qhull_parsing.hpp>
#include <scomplex/path_snapper.hpp>
#include <scomplex/chain-calc.hpp>

#include "boost/graph/graph_traits.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <tuple>

typedef std::vector<double> point_t;
typedef std::vector<size_t> cell_t;

int f(point_t a) { return a.at(0) < .5 ? 0 : 1; }

int main() {
    std::vector<point_t> points_v;
    std::vector<cell_t> cells_v;

    std::string filename;
    std::cout << "insert file name: ";
    std::cin >> filename;
    std::tie(points_v, cells_v) = parse_qhull_file(filename);

    simplicial::SimplicialComplex rsc(points_v, cells_v);
    auto qsc(rsc.quotient(f));
    snap::path_snapper snapper(rsc);

    auto path = snapper.snap_path({{-1, -1}, {1, 0}, {1, 1}});
    std::copy(path.begin(), path.end(),
              std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << '\n';

    bounding_chain chain_calc(points_v, cells_v);

    path = snapper.snap_path({{-1, -1}, {1, 0}, {1, 1}, {-1, -1}});
    std::copy(path.begin(), path.end(),
              std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << '\n';
    vector_t path_v(
        snapper.get_chain_vector({{-1, -1}, {1, 0}, {1, 1}, {-1, -1}}));

    // try {
    chain_t chain(2, path_v);
    chain = chain_calc.get_bounding_chain(chain);
    std::cout << path_v << '\n';
    vector_t b_chain = round_vec(std::get<1>(chain));
    std::cout << b_chain;
    // std::cout << chain_r << '\n';
    std::cout << round_vec(rsc.get_boundary(1) * std::get<1>(chain)) << '\n';
    std::cout << round_vec(rsc.get_boundary(1) * b_chain) << '\n';

    //} catch (NonZeroChain) {
    std::cout << "not a zero chain";
    //}

    return 0;
}
