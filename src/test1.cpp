/**
 **/
#include <iostream>
#include <iterator>

#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/qhull_parsing.hpp>
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

    auto s_comp=std::shared_ptr<gsimp::simplicial_complex>(new gsimp::simplicial_complex(points_v, cells_v));

    std::cout << "complex dimension: " << s_comp->dimension() << "\n";
    for (int i = 0; i <= s_comp->dimension(); ++i) {
        std::cout << "level " << i << " has " << s_comp->get_level_size(i)
                  << " cells\n";
    }

    gsimp::path_snapper p_snap(s_comp);
    std::vector<size_t> snapped = p_snap.snap_path_to_indices({{-2,-2},{-2,2},{2,2},{-2,-2}});
    std::copy(snapped.begin(),snapped.end(),std::ostream_iterator<size_t>(std::cout," "));

    return 0;
}
