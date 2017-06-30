/**
 * SIMPLICIAL COMPLEX TEST:
 *     takes as argument a path to a file containing the description of a
 *   triangulation as output by qhull (the last column of the point coordinates
 *   is ignored).
 *     constantly asks the user to input a number which will be interpreted as a
 *   1-cell index, and outputs the cofaces of the given cell, (together with the
 *   boundary index) and the faces of the cell.
 **/
#include <vector>
#include <iostream>
#include <iterator>

#include <scomplex/types.hpp>
#include <scomplex/qhull_parsing.hpp>
#include <scomplex/simplicial_complex.hpp>



int main(int argc, char* argv[]) {
    std::vector<gsimp::point_t> points_v;
    std::vector<gsimp::cell_t> cells_v;

    std::string filename;
    if (argc == 1) {
        std::cout << "insert file name: ";
        std::cin >> filename;
    } else
        filename = argv[1];
    std::tie(points_v, cells_v) = parse_qhull_file(filename);

    gsimp::simplicial_complex s_comp(points_v, cells_v);
    s_comp.calculate_hasse();


    std::cout << "done";

    return 0;
}
