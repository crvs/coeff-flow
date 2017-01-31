/**
 * @file Simplicial_Complex.h
 * @brief  utilities to create and manipulate geometric simplicial complexes
 * @author Joao Carvalho
 * @version 0.1
 * @date 2017-01-20
 *
 * Copyright (C)
 * 2017 - Joao Carvalho
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */
#ifndef QUOT_SIMPLICIAL_COMPLEX
#define QUOT_SIMPLICIAL_COMPLEX

#include <Eigen/Sparse>

#include <gudhi/Simplex_tree/Simplex_tree_siblings.h>
#include <gudhi/Simplex_tree.h>

#include <scomplex/base_utils.hpp>

#include <math.h>
#include <vector>
#include <set>

#include <typeinfo>

// point_t needs to have an empty constructor, and a size() function call,
// further it needs a begin() function call. Mostly it needs to behave like a
// std::vector

namespace simplicial {

typedef std::pair<int, std::vector<double>> chain_t;

template <typename point_t>
class SimplicialComplex {
   private:
    // the type of matrix to be used
    typedef typename Eigen::SparseMatrix<double> matrix_t;

    struct SimpleOptions : Gudhi::Simplex_tree_options_full_featured {
        // simplex tree options
        typedef int Vertex_handle;
    };
    typedef typename Gudhi::Simplex_tree<SimpleOptions> SimplexTree;
    typedef typename Gudhi::Simplex_tree<SimpleOptions>::Simplex_handle
        Simplex_handle;

    typedef std::vector<Simplex_handle *> level_t;
    typedef std::vector<level_t> levels_t;

   public:
    std::vector<std::vector<int>> get_level(int dimen) {
        std::vector<std::vector<int>> level;
        for (auto simp : simplices.complex_simplex_range()) {
            if (simplices.dimension(simp) == dimen) {
                std::vector<int> v_simp;
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

    SimplicialComplex<point_t>() { geometric_q = true; }

    /**
     * @brief  builds a simplicial complex object out of points and "triangles"
     *
     * @param points_a  points
     * @param tris      the top dimensional simplicial complex structure
     */
    SimplicialComplex<point_t>(std::list<point_t> points_a,
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
            // error (something is out-of-range)
            // calculate_matrices();
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

    SimplicialComplex<point_t> quotient(int f(point_t)) {
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

        auto quotient_sc = SimplicialComplex<point_t>(q_points, simp_list);
        quotient_sc.geometric_q = geometric;

        return quotient_sc;
    }

    /**
     * @brief  function deciding if a simplicial complex is geometric, note
     * that
     * quotient simplicial complexes are not geometric.
     */
    bool is_geometric() { return geometric_q; }

   private:
    bool geometric_q;
    std::vector<int> dimensions;

    std::vector<matrix_t> boundary_matrices;

    std::set<int> vertex_set(Simplex_handle s) {
        std::set<int> v_set;
        for (auto v : simplices.simplex_vertex_range(s)) {
            v_set.insert(s);
        }
        return v_set;
    }

    // need to do the boundary inclusion crap (you know what I mean)
    int boundary_index(Simplex_handle s_1, Simplex_handle s_2) {
        // s1 has to include into s2

        /*
        int dim1 = simplices.dimension(s_1);
        int dim2 = simplices.dimension(s_2);
        decltype(s_1) bdry_h;
        decltype(s_2) simp_h;
        if (dim2 == dim1 - 1) {
            simp_h = decltype(s_1)(s_1);
            bdry_h = decltype(s_2)(s_2);
        } else if (dim1 == dim2 - 1) {
            simp_h = decltype(s_2)(s_2);
            bdry_h = decltype(s_1)(s_1);
        } else {
            std::cerr << "CAN'T BE A BOUNDARY RELATION";
            throw;
        }

        bool found = false;
        int index = dim1;

        // <--- this is not working need to find the parent!!
        // simp_h = simp_h.parent();

        auto someth = simp_h->second.children()->parent();
        auto someth2 = someth.first;

        std::cout << "simp_h: " << typeid(simp_h).name() << std::endl;
        std::cout << "simp_h-fst: " << someth2 << std::endl;
        std::cout << "someth: " << typeid(someth).name() << std::endl;
        // std::cout << "someth-fst: " << someth2 << std::endl;

        while (not found) { // <-- this is not working yet
        if (simp_h == bdry_h) {found = true; index--;} else {index--;simp_h =
        simp_h->parent;bdry_h=bdry_h->parent;}
        }
        */

        auto it_1 = simplices.simplex_vertex_range(s_1).begin();
        auto end_1 = simplices.simplex_vertex_range(s_1).end();

        auto it_2 = simplices.simplex_vertex_range(s_2).begin();
        auto end_2 = simplices.simplex_vertex_range(s_2).end();

        int orient = 0;
        int unmatch = 0;

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

};      // namespace simplicial
#endif  // SIMPLICIAL_COMPLEX
