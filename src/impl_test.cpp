#include <iostream>

#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/qhull_parsing.hpp>

int my_char(gsimp::point_t point) { return 1; }

int main(int argc, char* argv[]) {
    std::vector<point_t> points_v;
    std::vector<cell_t> cells_v;

    std::string filename;
    if (argc == 1) {
        std::cout << "insert file name: ";
        std::cin >> filename;
    } else
        filename = argv[1];
    std::tie(points_v, cells_v) = parse_qhull_file(filename);

    gsimp::simplicial_complex s_comp(points_v, cells_v);

    gsimp::simplicial_complex q_comp = s_comp.quotient(my_char);

    std::cout << s_comp.dimension() << "\n";

    return 0;
}
