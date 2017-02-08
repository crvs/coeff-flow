#pragma once

#include <Eigen/Sparse>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Simplex_tree/Simplex_tree_siblings.h>

#include <scomplex/base_utils.hpp>

#include <set>
#include <vector>
#include <math.h>

typedef std::vector<double> point_t;
namespace simplicial {

class SimplicialComplex {
   private:
    // the type of matrix to be used
    typedef typename Eigen::SparseMatrix<double> matrix_t;
    typedef typename Eigen::SparseVector<double> chain_t;

    // simplex tree options
    struct SimpleOptions : Gudhi::Simplex_tree_options_full_featured {
        typedef size_t Vertex_handle;
    };
    typedef typename Gudhi::Simplex_tree<SimpleOptions> SimplexTree;
    typedef typename Gudhi::Simplex_tree<SimpleOptions>::Simplex_handle
        Simplex_handle;

    typedef std::vector<Simplex_handle *> level_t;
    typedef std::vector<level_t> levels_t;
    typedef std::vector<size_t> simplex_t;

   public:
    std::vector<simplex_t> get_level(int dimen) {
        std::vector<simplex_t> level;
        for (auto simp : simplices.complex_simplex_range()) {
            if (simplices.dimension(simp) == dimen) {
                simplex_t v_simp;
                for (auto v : simplices.simplex_vertex_range(simp)) {
                    v_simp.push_back(v);
                }
                level.push_back(v_simp);
            }
        }
        return level;
    }

    std::vector<point_t> points;
    SimplexTree simplices;

    int get_dimension(int level) { return dimensions.at(level); }

    SimplicialComplex() { geometric_q = true; }

    /**
     * @brief  builds a simplicial complex object out of points and "triangles"
     *
     * @param points_a  points
     * @param tris      the top dimensional simplicial complex structure
     */
    SimplicialComplex(std::list<point_t> points_a,
                      std::list<std::list<int>> tris) {
        geometric_q = true;
        points = std::vector<point_t>();
        for (auto pt : points_a) {
            points.push_back(pt);
        }
        for (auto s : tris) {
            // cant have simplices with repeated vertices so we dedupe the lists
            simplices.insert_simplex_and_subfaces(dedupe_list(s));
            int d = s.size() - 1;
            if (simplices.dimension() < d) {
                simplices.set_dimension(d);
            }
        }

        // we need to count i--simplices in order to get a key for them
        std::vector<int> count = std::vector<int>(simplices.dimension() + 1, 0);
        levels_t levels;
        levels = levels_t(simplices.dimension() + 1, level_t());

        auto simplex_range = simplices.complex_simplex_range();
        for (auto s : simplex_range) {
            int d = simplices.dimension(s);
            simplices.assign_key(s, count.at(d)++);
            levels.at(d).push_back(&s);
            auto v_range = simplices.simplex_vertex_range(s);
            auto sbd = simplices.boundary_simplex_range(s);
        }
        dimensions = count;
    }

    matrix_t get_boundary(int d) {
        if (boundary_matrices.size() == 0) {
            calculate_matrices();
        }
        if (d < boundary_matrices.size()) {
            return boundary_matrices.at(d);
        } else {
            return matrix_t(0, 0);
        }
    }

    SimplicialComplex quotient(int f(point_t)) {
        int n_points = 1;
        auto corresp = std::vector<int>();

        // making the correspondence between points in the original complex and
        // points in the quotient. point index 0 corresponds to the "virtual"
        // point (in case it exists).
        std::list<point_t> q_points;
        bool geometric = true;

        for (auto p : points) {
            if (f(p) != 0) {
                corresp.push_back(n_points++);
                q_points.push_back(p);
            } else {
                if (geometric) {
                    // need to realize that it is not geometric, (flip geometric
                    // and add the virtual point
                    auto pos = q_points.begin();
                    q_points.insert(pos, point_t());
                    geometric = false;
                }
                corresp.push_back(0);
            }
        }
        if (geometric) {
            for (int i = 0; i < corresp.size(); i++) {
                corresp.at(i) -= 1;
            }
        }

        // producing list of simplices
        std::list<std::list<int>> simp_list;

        for (auto s : simplices.complex_simplex_range()) {
            std::list<int> s_q;
            for (int v : simplices.simplex_vertex_range(s)) {
                s_q.push_back(corresp.at(v));
            }
            simp_list.push_back(dedupe_list(s_q));
        }

        auto quotient_sc = SimplicialComplex(q_points, simp_list);
        quotient_sc.geometric_q = geometric;

        return quotient_sc;
    }

    /**
     * @brief  function deciding if a simplicial complex is geometric, note
     * that
     * quotient simplicial complexes are not geometric.
     */
    bool is_geometric() { return geometric_q; }

    size_t get_simplex_index(simplex_t simp) {
        auto sh = simplices.find(simp);
        return simplices.key(sh);
    }

   private:
    bool geometric_q;
    std::vector<int> dimensions;

    std::vector<matrix_t> boundary_matrices;

    std::set<size_t> vertex_set(Simplex_handle s) {
        std::set<size_t> v_set;
        for (auto v : simplices.simplex_vertex_range(s)) {
            v_set.insert(v);
        }
        return v_set;
    }

    // calculate the index of s_1 in the boundary of s_2
    int boundary_index(Simplex_handle s_1, Simplex_handle s_2) {
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
                if (it_1 != end_1) {
                    it_1++;
                }
                it_2++;
            }
        }
        return pow(-1, orient);
    }

    /**
     * @brief  calculate the boundary matrices
     */
    void calculate_matrices() {
        boundary_matrices = std::vector<matrix_t>();
        for (int k = 0; k < simplices.dimension(); k++) {
            boundary_matrices.push_back(
                matrix_t(dimensions.at(k), dimensions.at(k + 1)));
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

    bool is_point(point_t pt) { return not(pt.size() == 0); }

};  // class SimplicialComplex
};  // namespace simplicial
