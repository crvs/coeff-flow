#include <vector>
#include <functional>
#include <Eigen/Sparse>

namespace gsimp {
// geometric types
typedef typename std::vector<double> point_t;
typedef typename std::vector<size_t> cell_t;

// linear algebra types
typedef typename Eigen::SparseMatrix<double> matrix_t;
typedef typename Eigen::SparseVector<double> vector_t;

// chains
typedef typename std::pair<int, vector_t> chain_t;

// characteristic functions
// typedef std::function<int(point_t)> char_fun_t;
};
