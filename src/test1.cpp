/**
 * PATH SNAPPING TEST:
 *     takes as argument a path to a file containing the description of a
 *   triangulation as output by qhull (the last column of the point coordinates
 *   is ignored).
 *     comes up with a path in this complex (hard coded) and snaps it to the
 *   complex, the only point of this is to test the libraries for errors.
 **/

#include <iostream>
#include <iterator>

#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/qhull_parsing.hpp>

// testing
#include <scomplex/path_snapper.hpp>

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

    auto s_comp = std::shared_ptr<gsimp::simplicial_complex>(
        new gsimp::simplicial_complex(points_v, cells_v));

    std::cout << "complex dimension: " << s_comp->dimension() << "\n";
    for (int i = 0; i <= s_comp->dimension(); ++i) {
        std::cout << "level " << i << " has " << s_comp->get_level_size(i)
                  << " cells\n";
    }

    gsimp::path_snapper p_snap(s_comp);
    std::vector<size_t> snapped =
        p_snap.snap_path_to_indices({{-2, -2}, {-2, 2}, {2, 2}, {-2, -2}});
    std::copy(snapped.begin(), snapped.end(),
              std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << '\n';

    auto chain = p_snap.index_sequence_to_chain(snapped);

    return 0;
}
