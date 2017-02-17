#pragma once

#include <Eigen/Sparse>
#include <functional>
#include <vector>


namespace gsimp {
// geometric types
typedef typename std::vector<double> point_t;
typedef typename std::vector<size_t> cell_t;

// linear algebra types
typedef typename Eigen::SparseMatrix<double> matrix_t;
typedef typename Eigen::SparseVector<double> vector_t;

// chains
typedef typename std::pair<int, vector_t> chain_t;
int& chain_dim(chain_t& p) { return std::get<0>(p); }
vector_t& chain_rep(chain_t& p) { return std::get<1>(p); }
double& chain_val(chain_t& p, size_t i) { return std::get<1>(p).coeffRef(i); }
}; // namespace gsimp
