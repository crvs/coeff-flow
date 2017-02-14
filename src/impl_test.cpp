#include <iostream>

#include <scomplex/simplicial_complex.hpp>
#include <scomplex/qhull_parsing.hpp>

int main() {
    std::vector<point_t> points_v;
    std::vector<cell_t> cells_v;

    std::string filename;
    std::cout << "insert file name: ";
    std::cin >> filename;
    std::tie(points_v, cells_v) = parse_qhull_file(filename);

    gsimp::simplicial_complex s_comp(points_v, cells_v);
    std::cout << s_comp.dimension();

    return 0;
}
