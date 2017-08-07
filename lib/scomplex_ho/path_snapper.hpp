#pragma once

#include <iostream>

#include <Eigen/Sparse>
#include <scomplex_ho/types.hpp>
#include <scomplex_ho/chains.hpp>
#include <scomplex_ho/simplicial_complex.hpp>
#include <scomplex_ho/nn_utils.hpp>
#include <scomplex_ho/graph_utils.hpp>

namespace gsimp {

class path_snapper {
   private:
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:
    // constructors and such
    path_snapper(std::vector<point_t>&, std::vector<cell_t>&);
    path_snapper(std::shared_ptr<simplicial_complex>);
    path_snapper(simplicial_complex&);
    ~path_snapper();
    path_snapper(path_snapper&);
    path_snapper& operator=(const path_snapper&);
    // the things we want to do
    std::vector<point_t> snap_path_to_points(std::vector<point_t>);
    std::vector<size_t> snap_path_to_indices(std::vector<point_t>);
    chain snap_path_to_dense_chain(std::vector<point_t>);
    chain snap_path_to_sparse_chain(std::vector<point_t>);
    chain snap_path_to_chain(std::vector<point_t>);
    // interconversion
    std::vector<point_t> index_sequence_to_point(std::vector<size_t>);
    std::vector<size_t> point_sequence_to_index(std::vector<point_t>);
    chain index_sequence_to_sparse_chain(std::vector<size_t>);
    chain index_sequence_to_dense_chain(std::vector<size_t>);
    chain point_sequence_to_chain(std::vector<point_t>);
    std::shared_ptr<simplicial_complex> get_underlying_complex();
};

struct path_snapper::impl {
    friend class clusterer;

    tree_t point_tree;     // defined in nn_utils.hpp
    graph_t vertex_graph;  // defined in graph_utils.hpp
    std::shared_ptr<simplicial_complex> s_comp;

    impl(std::shared_ptr<simplicial_complex> p_sc) {
        s_comp = p_sc;
        auto points = s_comp->get_points();
        make_tree(point_tree, points);
        vertex_graph = calculate_one_skelleton_graph(*s_comp);
    };

    impl(simplicial_complex& sc) {
        s_comp.reset(new simplicial_complex(sc));
        auto points = s_comp->get_points();
        make_tree(point_tree, points);
        vertex_graph = calculate_one_skelleton_graph(*s_comp);
    }

    impl(std::vector<point_t>& pts, std::vector<cell_t>& cells) {
        s_comp = std::shared_ptr<simplicial_complex>(
            new simplicial_complex(pts, cells));
        make_tree(point_tree, pts);
        vertex_graph = calculate_one_skelleton_graph(*s_comp);
    }

    ~impl(){};

    std::vector<size_t> snap_path(std::vector<point_t> path) {
        auto way_points = snap_points_to_indexes(point_tree, path);
        // std::cout << "processed path: \n";
        // for (size_t i = 0 ; i < way_points.size() ; ++i ) {
        //     for (auto c : path[i] ) std::cout << c << " ";
        //     auto pt = s_comp->get_point( way_points[i] );
        //     std::cout << " -> ";
        //     for (auto c : pt ) std::cout << c << " ";
        //     std::cout << "\n";
        // }
        auto snapped_path = complete_path(vertex_graph, way_points);
        return snapped_path;
    }

    std::vector<std::pair<cell_t, int>> index_pairs(std::vector<point_t> path) {
        auto vertex_path = snap_path(path);
        std::vector<std::pair<cell_t, int>> pair_seq;
        for (auto it = vertex_path.begin(); std::next(it) != vertex_path.end();
             ++it) {
            //
            //
            size_t src(*it), trg(*(std::next(it)));            //
            if (src < trg) {                                   //
                pair_seq.push_back(                            //
                    std::make_pair<std::vector<size_t>, int>(  //
                        std::vector<size_t>{src, trg}, 1));    //
            } else {                                           //
                pair_seq.push_back(                            //
                    std::make_pair<std::vector<size_t>, int>(  //
                        std::vector<size_t>{trg, src}, -1));   //
            }                                                  //
        }
        return pair_seq;
    }
};

path_snapper::path_snapper(std::shared_ptr<simplicial_complex> sc) {
    p_impl = std::shared_ptr<impl>(new impl(sc));
};

path_snapper::path_snapper(simplicial_complex& sc) {
    p_impl = std::shared_ptr<impl>(new impl(sc));
}

path_snapper::path_snapper(std::vector<point_t>& pts,
                           std::vector<cell_t>& cells) {
    p_impl = std::shared_ptr<impl>(new impl(pts, cells));
}

path_snapper::~path_snapper() {}

path_snapper::path_snapper(path_snapper& other) { p_impl = other.p_impl; }

path_snapper& path_snapper::operator=(const path_snapper& other) {
    p_impl = other.p_impl;
    return *this;
}

std::vector<point_t> path_snapper::snap_path_to_points(
    std::vector<point_t> path) {
    auto index_path = p_impl->snap_path(path);
    std::vector<point_t> point_path;
    for (size_t p : index_path)
        point_path.push_back(p_impl->s_comp->get_points().at(p));
    return point_path;
}

std::vector<size_t> path_snapper::snap_path_to_indices(
    std::vector<point_t> path) {
    return p_impl->snap_path(path);
}

chain path_snapper::snap_path_to_dense_chain(std::vector<point_t> path) {
    auto pair_seq = p_impl->index_pairs(path);
    chain rep = p_impl->s_comp->new_dense_chain(1);
    for (auto pair : pair_seq) {
        auto cell = std::get<0>(pair);
        if (cell[0] != cell[1]) {
            size_t index = p_impl->s_comp->cell_to_index(std::get<0>(pair));
            rep[index] += std::get<1>(pair);
        }
    }
    return rep;
}

chain path_snapper::snap_path_to_sparse_chain(std::vector<point_t> path) {
    auto pair_seq = p_impl->index_pairs(path);
    chain rep = p_impl->s_comp->new_sparse_chain(1);
    for (auto pair : pair_seq) {
        auto cell = std::get<0>(pair);
        if (cell[0] != cell[1]) {
            size_t index = p_impl->s_comp->cell_to_index(std::get<0>(pair));
            rep[index] += std::get<1>(pair);
        }
    }
    return rep;
}

chain path_snapper::snap_path_to_chain(std::vector<point_t> path) {
    return snap_path_to_sparse_chain(path);
}

std::vector<point_t> path_snapper::index_sequence_to_point(
    std::vector<size_t> ind_path) {
    std::vector<point_t> pt_path;
    for (size_t ind : ind_path)
        pt_path.push_back(p_impl->s_comp->get_points().at(ind));
    return pt_path;
}

std::vector<size_t> path_snapper::point_sequence_to_index(
    std::vector<point_t> pt_path) {
    std::vector<size_t> ind_path =
        snap_points_to_indexes(p_impl->point_tree, pt_path);
    return ind_path;
}

chain path_snapper::index_sequence_to_dense_chain(std::vector<size_t> ind_path) {
    auto it = ind_path.begin();
    chain rep = p_impl->s_comp->new_dense_chain(1);
    for (; std::next(it) != ind_path.end(); ++it) {
        auto ind = p_impl->s_comp->cell_to_index({*it, *std::next(it)});
        rep[ind] += (*it < *std::next(it)) ? 1 : -1;
    }
    return rep;
}

chain path_snapper::index_sequence_to_sparse_chain(std::vector<size_t> ind_path) {
    auto it = ind_path.begin();
    chain rep = p_impl->s_comp->new_chain(1);
    for (; std::next(it) != ind_path.end(); ++it) {
        auto ind = p_impl->s_comp->cell_to_index({*it, *std::next(it)});
        rep[ind] += (*it < *std::next(it)) ? 1 : -1;
    }
    return rep;
}

chain path_snapper::point_sequence_to_chain(std::vector<point_t> pt_path) {
    return index_sequence_to_sparse_chain(point_sequence_to_index(pt_path));
}

std::shared_ptr<simplicial_complex> path_snapper::get_underlying_complex() {
    return p_impl->s_comp;
}

};
