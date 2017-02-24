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
#include <scomplex/drawing_utils.hpp>
#include <scomplex/utils.hpp>

// testing chain_calc
#include <scomplex/chain_calc.hpp>

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

    std::vector<point_t> path2 = {{-2, -2}, { 2, -2}, { 2,  2}, {-2, -2}};
    std::vector<point_t> path  = {{-2, -2}, {.3, .3}, {.2, -2}, {-2, -2}};

    // auxiliary structures
    gsimp::path_snapper p_snap(s_comp);
    gsimp::bounding_chain c_calc(s_comp);

    gsimp::draw_path(p_snap, drawing, path, 0.003);
    gsimp::draw_path(p_snap, drawing, path2, 0.003);

    gsimp::chain_t path_chain = p_snap.snap_path_to_chain(path);
    gsimp::add_to(p_snap.snap_path_to_chain(path2),path_chain);
    try {
        gsimp::chain_t bound_chain = c_calc.get_bounding_chain(path_chain);
        draw_2_chain(s_comp, drawing, bound_chain,
                     "#ff00ff' fill-opacity='0.8");
    } catch (const std::exception& e) {
        std::cout << "oh well, nothing I can do";
    }

    return 0;
}
