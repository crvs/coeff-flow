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

    for (auto cell : s_comp.get_level(2)) {
        std::vector<std::vector<double>> cell_pt;
        for (auto i : cell) cell_pt.push_back(points_v[i]);
        drawing.draw_triangle(cell_pt, "#00ff00");
    }

    for (auto cell : s_comp.get_level(1)) {
        std::vector<std::vector<double>> cell_pt;
        for (auto i : cell) cell_pt.push_back(points_v[i]);
        drawing.draw_edge(cell_pt[0], cell_pt[1], "#000000", 0.002);
    }

    gsimp::path_snapper p_snap(s_comp);
    std::vector<size_t> snapped =
        p_snap.snap_path_to_indices({{-2, -2}, {-2, 2}, {2, 2}, {-2, -2}});
    auto it = snapped.begin();
    auto trail = it++;
    for (; it != snapped.end(); ++it, ++trail) {
        auto init = s_comp.get_points()[*trail];
        auto dest = s_comp.get_points()[*it];
        drawing.draw_edge(init, dest, "#ff0000", 0.003);
    }

    return 0;
}
