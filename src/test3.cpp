/**
 * PATH SNAPPING TEST:
 *     takes as argument a path to a file containing the description of a
 *   triangulation as output by qhull (the last column of the point coordinates
 *   is ignored).
 *     outputs the complex to an svg file called my_complex
 **/

#include <iostream>
#include <iterator>

#include <scomplex/types.hpp>
#include <scomplex/qhull_parsing.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/path_snapper.hpp>
#include <scomplex/utils.hpp>

// testing
#include <scomplex/drawing_utils.hpp>

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

    gsimp::drawing drawing("my_complex", 500, 500);

    gsimp::draw_complex(s_comp, drawing, 0.003, "#00ff00", "#000000");

    { // scope for p_snap
        gsimp::path_snapper p_snap(s_comp);
        gsimp::draw_path(p_snap, drawing, {{-2, -2}, {2, -2}, {2, 2}, {-2, -2}},
                         0.003);
    }

    return 0;
}
