/**
 * SIMPLICIAL COMPLEX TEST:
 *     takes as argument a path to a file containing the description of a
 *   triangulation as output by qhull (the last column of the point coordinates
 *   is ignored).
 *     constantly asks the user to input a number which will be interpreted as a
 *   1-cell index, and outputs the cofaces of the given cell, (together with the
 *   boundary index) and the faces of the cell.
 **/
#include <iostream>
#include <iterator>

#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/qhull_parsing.hpp>

int my_char(gsimp::point_t point) { return 1; }

int& do_this(int& a) { return a; }

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

    std::cout << "complex dimension: " << s_comp.dimension() << "\n";
    for (int i = 0; i <= s_comp.dimension(); ++i) {
        std::cout << "level " << i << " has " << s_comp.get_level_size(i)
                  << " cells\n";
    }

    size_t i;
    std::cout << "get 1-cell number: ";
    while (std::cin >> i) {
        std::vector<std::pair<int, size_t>> my_cofaces =
            s_comp.get_cof_and_ind_index(1, i);
        for (auto s : my_cofaces) {
            std::cout << std::get<0>(s) << ", ";
            std::cout << std::get<1>(s) << " ;  ";
        }
        std::cout << '\n';
        std::vector<size_t> faces = s_comp.cell_boundary_index(1, i);
        std::copy(faces.begin(), faces.end(),
                  std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << "\nget 1-cell number: ";
    }

    return 0;
}
