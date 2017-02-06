#pragma once

namespace snap {

#include <scomplex/nn_utils.hpp>
#include <scomplex/graph_utils.hpp>
#include <scomplex/Simplicial_Complex.h>

class path_snapper {
   private:
    tree_t point_tree;   // defined in nn_utils.hpp
    Graph vertex_graph;  // defined in graph_utils.hpp

   public:
    path_snapper(simplicial::SimplicialComplex<point_t>& sc) {
        make_tree(point_tree, sc.points);
        vertex_graph = calculate_one_skelleton_graph(sc);
    }

    std::vector<size_t> snap_path(std::vector<point_t> path) {
        std::vector<size_t> way_points =
            snap_points_to_indexes(point_tree, path);
        std::copy(way_points.begin(), way_points.end(),
                  std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << '\n';

        std::vector<size_t> snapped_path =
            complete_path(vertex_graph, way_points);
        return snapped_path;
    };
};
};
