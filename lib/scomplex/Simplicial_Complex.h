#ifndef SIMPLICIAL_COMPLEX
#define SIMPLICIAL_COMPLEX

#include <string.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gudhi/Simplex_tree/Simplex_tree_siblings.h>
#include <gudhi/Simplex_tree.h>

#include <math.h>
#include <vector>
#include <set>

// point_t needs to have an empty constructor, and a size() function call,
// further it needs a begin() function call. Mostly it needs to behave like a
// std::vector

template <typename point_t>
class SimplicialComplex {
   private:
    typedef typename boost::numeric::ublas::compressed_matrix<int> matrix_t;

    struct SimpleOptions : Gudhi::Simplex_tree_options_full_featured {
        // simplex tree options
        typedef int Vertex_handle;
    };
    typedef typename Gudhi::Simplex_tree<SimpleOptions> SimplexTree;
    typedef typename Gudhi::Simplex_tree<SimpleOptions>::Simplex_handle
        Simplex_handle;

   public:
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

    std::list<point_t> points;
    SimplexTree simplices;

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
        points = points_a;
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

        auto simplex_range = simplices.complex_simplex_range();
        for (auto s : simplex_range) {
            int d = simplices.dimension(s);
            simplices.assign_key(s, count.at(d)++);
            auto v_range = simplices.simplex_vertex_range(s);
            auto sbd = simplices.boundary_simplex_range(s);
        }
        dimensions = count;
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
                    auto pos = points.begin();
                    points.insert(pos, point_t());
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

    std::string print_level(int level) {
        std::string out;
        return out;
    }

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
    int boundary_inclusion_orientation(Simplex_handle s_1, Simplex_handle s_2) {
        // s1 includes into s2
        auto it_1 = simplices.simplex_vertex_range(s_1).begin();
        auto it_2 = simplices.simplex_vertex_range(s_2).begin();
        int orient = 0;
        for (int p = 0; p <= simplices.dimension(s_1); p++) {
            if (it_1 != it_2) {
                it_2++;
                orient = p;
            } else {
                it_1++;
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
                boundary_matrices.at(k)(i, j) =
                    boundary_inclusion_orientation(bs, s);
            }
        }
    }

    int quick_dim(Simplex_handle s) const {
        int s_dim = 0;
        for (auto v : simplices.simplex_vertex_range(s)) {
            s_dim++;
        }
    }

    bool is_point(point_t pt) { return not(pt.size() == 0); }
    static std::list<int> dedupe_list(std::list<int> list) {
        // O(n log(n))
        std::set<int> no_reps;
        std::list<int> no_reps_list;
        for (int el : list) {
            no_reps.insert(el);
        }
        for (int el : no_reps) {
            no_reps_list.push_back(el);
        }
        return no_reps_list;
    }
};

#endif
