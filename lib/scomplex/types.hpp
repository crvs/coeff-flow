#pragma once

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace gsimp {
// geometric types
typedef std::vector<double> point_t;
typedef std::vector<size_t> cell_t;

// linear algebra types
typedef typename Eigen::SparseMatrix<double> matrix_t;
typedef typename Eigen::SparseVector<double> vector_t;

// chains
typedef typename std::pair<int, vector_t> chain_t;

int& chain_dim(chain_t& p) { return std::get<0>(p); }
vector_t& chain_rep(chain_t& p) { return std::get<1>(p); }
double& chain_val(chain_t& p, size_t i) { return std::get<1>(p).coeffRef(i); }

chain_t add(chain_t chain1, chain_t chain2) {
    if (std::get<0>(chain1) == std::get<0>(chain2))
        return chain_t(std::get<0>(chain1),
                       std::get<1>(chain1) + std::get<1>(chain2));
    std::cout << "chains must have the same dimension to be added (+)";
    throw std::exception();
}

void add_to(chain_t chain1, chain_t chain2) {
    if (std::get<0>(chain1) == std::get<0>(chain2))
        std::get<1>(chain1) += std::get<1>(chain2);
    else {
        std::cout << "chains must have the same dimension to be added (+=)";
        throw std::exception();
    }
}

chain_t subtract(chain_t chain1, chain_t chain2) {
    if (std::get<0>(chain1) == std::get<0>(chain2))
        return chain_t(std::get<0>(chain1),
                       std::get<1>(chain1) - std::get<1>(chain2));
    std::cout << "chains must have the same dimension to be added (-)";
    throw std::exception();
}

void subtract_to(chain_t chain1, chain_t chain2) {
    if (std::get<0>(chain1) == std::get<0>(chain2))
        std::get<1>(chain1) -= std::get<1>(chain2);
    else {
        std::cout << "chains must have the same dimension to be added (-=)";
        throw std::exception();
    }
}

chain_t prod(double coef, chain_t chain) {
    return chain_t(std::get<0>(chain), std::get<1>(chain) * coef);
}

void prod_to(double coef, chain_t chain) { std::get<1>(chain) *= coef; }

Eigen::VectorXd point_to_eigen(point_t pt) {
    Eigen::VectorXd v0(pt.size());
    for (int j = 0; j < pt.size(); ++j) v0(j) = pt[j];
    return v0;
}
}  // namespace gsimp
