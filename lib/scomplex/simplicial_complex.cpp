#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>
#include <scomplex/base_utils.hpp>

#include <iterator>  // for debuging purposes

#include <gudhi/Simplex_tree.h>
#include <Eigen/Sparse>

#include <tuple>
#include <memory>  // smart pointers
#include <functional>
#include <cmath>

namespace gsimp {

struct simplicial_complex::impl {
    // auxiliary types
    struct SimpleOptions : Gudhi::Simplex_tree_options_full_featured {
        typedef size_t Vertex_handle;
    };
    typedef Gudhi::Simplex_tree<SimpleOptions> simp_tree;
    typedef Gudhi::Simplex_tree<SimpleOptions>::Simplex_handle simp_handle;
    typedef std::vector<simp_handle*> level_t;
    typedef std::vector<std::unique_ptr<level_t>> levels_t;

    // member variables
    std::vector<point_t> points;
    simp_tree simplices;
    bool quotient_q;
    std::vector<matrix_t> boundary_matrices;
    levels_t levels;

    impl() : quotient_q(false){};

    impl(std::vector<point_t>& arg_points, std::vector<cell_t>& arg_tris)
        : points(arg_points), quotient_q(false) {
        // create the simplex tree
        for (auto tri : arg_tris) {
            simplices.insert_simplex_and_subfaces(dedupe_vec(tri));
            int d = tri.size() - 1;
            if (simplices.dimension() < d) simplices.set_dimension(d);
        }

        // assign a key to each simplex in each level
        std::vector<size_t> count(simplices.dimension() + 1, 0);
        for (int i = 0; i < simplices.dimension() + 1; ++i) {
            level_t* level = new level_t();
            levels.push_back(std::unique_ptr<level_t>(level));
        }

        auto simplex_range = simplices.complex_simplex_range();
        for (auto s : simplex_range) {
            int d = simplices.dimension(s);
            simplices.assign_key(s, count.at(d)++);
            levels.at(d)->push_back(new simp_handle(s));
        }
    }

    ~impl() {
        // get rid of the levels so they won't dangle
        for (int i = levels.size(); i <= 0; --i) levels.at(i).reset();
    };

    size_t get_level_size(int level) { return levels.at(level)->size(); }

    // calculate the index of s_1 in the boundary of s_2
    int boundary_index(simp_handle s_1, simp_handle s_2) {
        auto it_1 = simplices.simplex_vertex_range(s_1).begin();
        auto end_1 = simplices.simplex_vertex_range(s_1).end();
        //
        auto it_2 = simplices.simplex_vertex_range(s_2).begin();
        //
        int orient = 0;
        int unmatch = 0;
        //
        for (int p = 0; p <= simplices.dimension(s_2); p++) {
            if (*it_1 != *it_2) {
                orient = p++;
                unmatch++;
            } else {
                if (it_1 != end_1) it_1++;
                it_2++;
            }
        }
        return pow(-1, orient);
    }

    void calculate_matrices() {
        boundary_matrices = std::vector<matrix_t>();
        for (int k = 0; k < simplices.dimension(); k++) {
            boundary_matrices.push_back(
                matrix_t(get_level_size(k), get_level_size(k + 1)));
        }
        for (auto s : simplices.complex_simplex_range()) {
            int j = simplices.key(s);
            for (auto bs : simplices.boundary_simplex_range(s)) {
                int i = simplices.key(bs);
                int k = simplices.dimension(bs);
                boundary_matrices.at(k).coeffRef(i, j) = boundary_index(bs, s);
            }
        }
    }

    std::vector<cell_t> get_level(int level) {
        std::vector<cell_t> level_cells;
        for (auto simp : *levels.at(level)) {
            cell_t v_simp;
            for (auto v : simplices.simplex_vertex_range(*simp)) {
                v_simp.push_back(v);
            }
            level_cells.push_back(v_simp);
        }
        return level_cells;
    }

    simp_handle index_to_handle(int d, size_t tau) {
        return *(levels.at(d)->at(tau));
    }

    size_t handle_to_index(simp_handle tau) { return simplices.key(tau); }

    simp_handle cell_to_handle(cell_t tau) {
        auto sh = simplices.find(tau);
        return sh;
    }

    cell_t handle_to_cell(simp_handle tau) {
        cell_t cell;
        for (auto v : simplices.simplex_vertex_range(tau)) cell.push_back(v);
        return cell;
    }
};  // struct impl

std::vector<std::pair<int, cell_t>> simplicial_complex::get_bdry_and_ind(
    cell_t cell) {
    impl::simp_handle simp = p_impl->cell_to_handle(cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<std::pair<int, cell_t>> boundary_and_indices;
    for (auto face : c_boundary)
        boundary_and_indices.push_back(              //
            std::make_pair<int, cell_t>(             //
                p_impl->boundary_index(face, simp),  //
                p_impl->handle_to_cell(face)));      //
    return boundary_and_indices;
};

std::vector<std::pair<int, size_t>> simplicial_complex::get_bdry_and_ind_index(
    int d, size_t cell) {
    impl::simp_handle simp = p_impl->index_to_handle(d, cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<std::pair<int, size_t>> boundary_and_indices;
    for (auto face : c_boundary)
        boundary_and_indices.push_back(              //
            std::make_pair<int, size_t>(             //
                p_impl->boundary_index(face, simp),  //
                p_impl->handle_to_index(face)));     //
    return boundary_and_indices;
};

std::vector<size_t> simplicial_complex::cell_boundary_index(int d,
                                                            size_t cell) {
    impl::simp_handle simp = p_impl->index_to_handle(d, cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<size_t> s_boundary;
    for (auto c : c_boundary) s_boundary.push_back(p_impl->handle_to_index(c));
    return s_boundary;
}

std::vector<cell_t> simplicial_complex::cell_boundary(cell_t cell) {
    impl::simp_handle simp = p_impl->cell_to_handle(cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<cell_t> s_boundary;
    for (auto c : c_boundary) s_boundary.push_back(p_impl->handle_to_cell(c));
    return s_boundary;
}

int simplicial_complex::boundary_inclusion_index(cell_t c1, cell_t c2) {
    return p_impl->boundary_index(p_impl->cell_to_handle(c1),   //
                                  p_impl->cell_to_handle(c2));  //
};
int simplicial_complex::boundary_inclusion_index(int d1, size_t s1,    //
                                                 int d2, size_t s2) {  //
    cell_t c1 = index_to_cell(d1, s1);
    cell_t c2 = index_to_cell(d2, s2);
    return p_impl->boundary_index(p_impl->cell_to_handle(c1),
                                  p_impl->cell_to_handle(c2));
};

std::vector<std::pair<int, cell_t>> simplicial_complex::get_cof_and_ind(
    cell_t cell) {
    std::vector<std::pair<int, cell_t>> c_cofaces;
    for (auto face : get_cofaces(cell))
        c_cofaces.push_back(                                   //
            std::make_pair(                                    //
                boundary_inclusion_index(cell, face), face));  //
    return c_cofaces;
}

std::vector<std::pair<int, size_t>> simplicial_complex::get_cof_and_ind_index(
    int d, size_t c) {
    std::vector<std::pair<int, size_t>> c_cofaces;
    for (auto face : get_cofaces_index(d, c))
        c_cofaces.push_back(                                          //
            std::make_pair(                                           //
                boundary_inclusion_index(d, c, d + 1, face), face));  //
    return c_cofaces;
}

int simplicial_complex::get_level_size(int level) {
    return p_impl->get_level_size(level);
}

simplicial_complex::simplicial_complex(std::vector<point_t>& arg_points,
                                       std::vector<cell_t>& arg_tris) {
    p_impl = std::shared_ptr<impl>(new impl(arg_points, arg_tris));
}

simplicial_complex::simplicial_complex(const simplicial_complex& other) {
    p_impl = other.p_impl;
}

simplicial_complex& simplicial_complex::operator=(
    const simplicial_complex& other) {
    p_impl = other.p_impl;
    return *this;
}

simplicial_complex::~simplicial_complex() {}

std::vector<point_t> simplicial_complex::get_points() { return p_impl->points; }

std::vector<cell_t> simplicial_complex::get_level(int level) {
    return p_impl->get_level(level);
}

matrix_t simplicial_complex::get_boundary_matrix(int d) {
    // uninstantiated boundary matrices
    if (p_impl->boundary_matrices.size() == 0) p_impl->calculate_matrices();

    // now they have to be instantiated, get them
    if (0 <= d && d < p_impl->boundary_matrices.size())
        return p_impl->boundary_matrices.at(d);
    else
        throw NoChain();
}

int simplicial_complex::dimension() { return p_impl->simplices.dimension(); }

cell_t simplicial_complex::index_to_cell(int d, size_t ind) {
    auto sh = p_impl->levels.at(d)->at(ind);
    return p_impl->handle_to_cell(*sh);
}
size_t simplicial_complex::cell_to_index(cell_t simp) {
    auto sh = p_impl->simplices.find(simp);
    return p_impl->simplices.key(sh);
}

bool simplicial_complex::is_quotient() { return p_impl->quotient_q; }

std::vector<size_t> simplicial_complex::get_cofaces_index(int d, size_t face) {
    std::vector<size_t> s_cofaces;
    impl::simp_handle face_h = *(p_impl->levels.at(d)->at(face));
    cell_t my_cell = p_impl->handle_to_cell(face_h);
    // codimension 1 faces
    auto range = p_impl->simplices.cofaces_simplex_range(face_h, 1);
    for (auto tau : range) s_cofaces.push_back(p_impl->simplices.key(tau));
    return s_cofaces;
}

std::vector<cell_t> simplicial_complex::get_cofaces(cell_t face) {
    std::vector<cell_t> s_cofaces;
    auto face_h = p_impl->cell_to_handle(face);
    // codimension 1 faces
    auto range = p_impl->simplices.cofaces_simplex_range(face_h, 1);
    for (auto tau : range) s_cofaces.push_back(p_impl->handle_to_cell(tau));
    return s_cofaces;
}

// TODO: get top dimensional cells from the simplicial complex to add to the
// simplex tree (instead of all the simplices).
simplicial_complex simplicial_complex::quotient(int char_fun(point_t)) {
    int n_points = 1;
    auto corresp = std::vector<size_t>();

    // making the correspondence between points in the original complex
    // and points in the quotient. point index 0 corresponds to the
    // "virtual" point (in case it exists).
    std::vector<point_t> points;
    bool quotient_q = p_impl->quotient_q;

    for (auto p : p_impl->points) {
        if (char_fun(p) != 0) {
            corresp.push_back(n_points++);
            points.push_back(p);
        } else {
            if (not quotient_q) {
                // need to realize that it is not geometric, (flip
                // quotient_q and add the virtual point
                points.insert(points.begin(), point_t());
                quotient_q = true;
            }
            corresp.push_back(0);
        }
    }

    if (not quotient_q)
        for (int i = 0; i < corresp.size(); i++) corresp.at(i) -= 1;

    // producing list of simplices
    std::vector<cell_t> simp_list;

    for (auto s : p_impl->simplices.complex_simplex_range()) {
        cell_t s_q;
        for (int v : p_impl->simplices.simplex_vertex_range(s))
            s_q.push_back(corresp.at(v));
        simp_list.push_back(s_q);
    }

    simplicial_complex quotient_sc(points, simp_list);
    quotient_sc.p_impl->quotient_q = quotient_q;

    return quotient_sc;
}

};  // namespace gsimp

// DEPRECATED

/*
 *std::set<size_t> vertex_set(simp_handle s) {
 *    std::set<size_t> v_set;
 *    for (auto v : simplices.simplex_vertex_range(s)) {
 *        v_set.insert(v);
 *    }
 *    return v_set;
 *}
 */

/* MAY BE INCLUDED AGAIN WHEN WE START DOING MORE WITH QUOTIENTS
 *
 * bool is_point(point_t pt) {
 *     return (not quotient_q ? true : not(pt.size() == 0));
 * }
 *
 */
