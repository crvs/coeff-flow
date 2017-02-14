#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>
#include <scomplex/base_utils.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Simplex_tree/Simplex_tree_siblings.h>
#include <Eigen/Sparse>
#include <memory>
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
        : quotient_q(false), points(arg_points) {
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
            levels.at(d)->push_back(&s);
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
        auto end_2 = simplices.simplex_vertex_range(s_2).end();
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
};  // struct impl

int simplicial_complex::get_level_size(int level) {
    return p_impl->get_level_size(level);
}

simplicial_complex::simplicial_complex(std::vector<point_t>& arg_points,
                                       std::vector<cell_t>& arg_tris) {
    p_impl = std::shared_ptr<impl>(new impl(arg_points, arg_tris));
}

simplicial_complex::~simplicial_complex() {}

std::vector<point_t> simplicial_complex::get_points() { return p_impl->points; }

std::vector<cell_t> simplicial_complex::get_level(int level) {
    return p_impl->get_level(level);
}

matrix_t simplicial_complex::get_boundary_matrix(int d) {
    // uninstantiated boundary matrices
    if (p_impl->boundary_matrices.size() == 0) {
        p_impl->calculate_matrices();
    }

    // now they have to be instantiated, get them
    if (0 <= d && d < p_impl->boundary_matrices.size()) {
        return p_impl->boundary_matrices.at(d);
    } else {
        throw NoChain();
    }
}

int simplicial_complex::dimension() { return p_impl->simplices.dimension(); }

size_t simplicial_complex::get_index_of_simplex(cell_t simp) {
    auto sh = p_impl->simplices.find(simp);
    return p_impl->simplices.key(sh);
}

bool simplicial_complex::is_quotient() { return p_impl->quotient_q; }

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
                quotient_q = false;
            }
            corresp.push_back(0);
        }
    }

    if (not quotient_q) {
        for (int i = 0; i < corresp.size(); i++) {
            corresp.at(i) -= 1;
        }
    }

    // producing list of simplices
    std::vector<cell_t> simp_list;

    for (auto s : p_impl->simplices.complex_simplex_range()) {
        cell_t s_q;
        for (int v : p_impl->simplices.simplex_vertex_range(s))
            s_q.push_back(corresp.at(v));
        simp_list.push_back(s_q);
    }

    auto quotient_sc = simplicial_complex(points, simp_list);
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
