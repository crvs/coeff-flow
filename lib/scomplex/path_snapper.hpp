#pragma once

namespace snap {

#include <scomplex/nn_utils.hpp>
#include <scomplex/graph_utils.hpp>
#include <scomplex/Simplicial_Complex.h>
#include <Eigen/Sparse>

typedef std::vector<std::pair<std::vector<size_t>, int>> geo_chain_t;
typedef Eigen::SparseVector<double> vector_t;

class path_snapper {
   private:
    tree_t point_tree;   // defined in nn_utils.hpp
    Graph vertex_graph;  // defined in graph_utils.hpp
    simplicial::SimplicialComplex<point_t>* s_comp;

   public:
    path_snapper(simplicial::SimplicialComplex<point_t>& sc) {
        make_tree(point_tree, sc.points);
        vertex_graph = calculate_one_skelleton_graph(sc);
        s_comp = &sc;
    }

    std::vector<size_t> snap_path(std::vector<point_t> path) {
        std::vector<size_t> way_points =
            snap_points_to_indexes(point_tree, path);
        std::vector<size_t> snapped_path =
            complete_path(vertex_graph, way_points);
        return snapped_path;
    }

    geo_chain_t get_chain_rep(std::vector<point_t> path) {
        auto vertex_path = snap_path(path);
        geo_chain_t chain_rep;
        for (                                    //
            auto it = vertex_path.begin();       //
            std::next(it) != vertex_path.end();  //
            ++it) {
            size_t s(*it), t(*(std::next(it)));
            if (s < t) {
                chain_rep.push_back(  //
                    std::make_pair<std::vector<size_t>, int>(
                        std::vector<size_t>{s, t}, 1));
            } else {
                chain_rep.push_back(  //
                    std::make_pair<std::vector<size_t>, int>(
                        std::vector<size_t>{s, t}, -1));
            }
        }
    }

    vector_t get_chain_vector(std::vector<point_t> path) {
        auto chain_rep = get_chain_rep(path);
        vector_t vector_rep(s_comp->get_dimension(1));
        for (auto chain_el : chain_rep) {
            simplex_t simp = std::get<0>(chain_el);
            int val = std::get<1>(chain_el);
            size_t index = s_comp->get_simplex_index(simp);
            vector_rep.coeffRef(index) = val;
        }
    }
};
};
