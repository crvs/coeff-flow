#pragma once

#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <cmath>
#include <exception>
#include <tuple>
#include <vector>

#include <memory>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <fstream>

namespace gsimp {

class non_zero_chain : public std::exception {};
class cant_draw_dimension : public std::exception {};

class bounding_chain {
    std::shared_ptr<simplicial_complex> s_comp;
    std::vector<std::unique_ptr<matrix_t>> boundary_matrices;
    void populate_matrices();

   public:
    bounding_chain(simplicial_complex& sc);
    bounding_chain(std::shared_ptr<simplicial_complex> sc);
    bounding_chain(std::vector<point_t>& points, std::vector<cell_t>& tris);
    ~bounding_chain();
    chain_t get_bounding_chain(chain_t&);
    // void draw_chain(std::string, chain_t);
};

//---------------------------------------
// implementation

bounding_chain::~bounding_chain() {}

bounding_chain::bounding_chain(std::shared_ptr<simplicial_complex> sc) {
    s_comp = sc;
    populate_matrices();
}

bounding_chain::bounding_chain(simplicial_complex& sc) {
    s_comp.reset(new simplicial_complex(sc));
    populate_matrices();
}

bounding_chain::bounding_chain(std::vector<point_t>& points,
                               std::vector<cell_t>& tris) {
    s_comp.reset(new simplicial_complex(points, tris));
    populate_matrices();
}

void bounding_chain::populate_matrices() {
    for (int d = 0; d < s_comp->dimension(); ++d) {
        auto level_matrix = s_comp->get_boundary_matrix(d);
        boundary_matrices.push_back(
            std::unique_ptr<matrix_t>(new matrix_t(level_matrix)));
    }
}

// round vectors for comparison
vector_t round_vec(vector_t vec) {
    vector_t rounded_vec(vec.rows());
    for (int i = 0; i < vec.rows(); ++i) {
        int coef = round(vec.coeffRef(i));
        if (coef != 0) rounded_vec.coeffRef(i) = coef;
    }
    return rounded_vec;
}

bool equals(vector_t vec1, vector_t vec2) {
    for (int i = 0; i < vec1.rows(); ++i) {
        if (vec1.coeffRef(i) != vec2.coeffRef(i)) {
            std::cout << i << ": " << vec1.coeffRef(i) << " "
                      << vec2.coeffRef(i) << '\n';
            return false;
        }
    }
    return true;
}

chain_t bounding_chain::get_bounding_chain(chain_t& chain) {
    int chain_d;
    vector_t chain_v;
    std::tie<int, vector_t>(chain_d, chain_v) = chain;

    if (chain_d >= s_comp->dimension()) throw non_zero_chain();

    Eigen::LeastSquaresConjugateGradient<matrix_t> lscg;
    lscg.compute(*(boundary_matrices.at(chain_d)));

    vector_t bound_chain(round_vec(lscg.solve(chain_v)));
    vector_t result(round_vec(*(boundary_matrices.at(chain_d)) * bound_chain));

    if (equals(result, chain_v))
        return chain_t(chain_d - 1, bound_chain);
    else
        throw non_zero_chain();
}

};
