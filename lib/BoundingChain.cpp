#include "BoundingChain.hpp"

#include "chains.hpp"
#include "types.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

/*
 * This is the implementation file to the header "BoundingChain.hpp" it
 * provides the facilities to calculate bounding chains for a given chain in on
 * a "SimplicialComplex" object.
 * */

namespace Gsimp {

BoundingChain::~BoundingChain() {} // all used pointer are smart

BoundingChain::BoundingChain(std::shared_ptr< SimpComplex > sc) : s_comp{sc} {
    populate_matrices();
}

BoundingChain::BoundingChain(SimpComplex& sc) : s_comp{new SimpComplex(sc)} {
    populate_matrices();
}

BoundingChain::BoundingChain(std::vector< Point >& points,
                             std::vector< Cell >& tris)
    : s_comp{new SimpComplex(points, tris)} {
    populate_matrices();
}

void BoundingChain::populate_matrices() {
    for (int d = 0; d < s_comp->dimension(); ++d) {
        auto level_matrix = s_comp->get_boundary_matrix(d);
        boundary_matrices.push_back(
            std::unique_ptr< Matrix >(new Matrix(level_matrix)));
    }
}

// round vectors for comparison
Vector round_vec(Vector vec) {
    Vector rounded_vec(vec.rows());
    for (int i = 0; i < vec.rows(); ++i) {
        int coef = round(vec.coeffRef(i));
        if (coef != 0) rounded_vec.coeffRef(i) = coef;
    }
    return rounded_vec;
}

bool equals(Vector vec1, Vector vec2) {
    for (int i = 0; i < vec1.rows(); ++i) {
        if (vec1.coeffRef(i) != vec2.coeffRef(i)) {
            return false;
        }
    }
    return true;
}

chain BoundingChain::getBoundingChain(chain& rep) {
    int chain_d = rep.dimension();
    Vector chain_v = rep.get_sparse();

    if (chain_d >= s_comp->dimension()) throw NonBoundary();

    Eigen::LeastSquaresConjugateGradient< Matrix > lscg;
    lscg.compute(*(boundary_matrices.at(chain_d)));

    Vector bound_chain(round_vec(lscg.solve(chain_v)));
    Vector result(round_vec(*(boundary_matrices.at(chain_d)) * bound_chain));

    if (equals(result, chain_v))
        return chain(chain_d + 1, bound_chain);
    else
        throw NonBoundary();
}

};
