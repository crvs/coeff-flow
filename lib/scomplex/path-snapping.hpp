#pragma once

namespace snap {

#include <scomplex/simplicial_complex.h>
#include <Eigen/Sparse>
#include <scomplex/graph_utils.hpp>
#include <scomplex/nn_utils.hpp>
#include <scomplex/types.hpp>

typedef Eigen::SparseVector<double> vector_t;

class path_snapper {
   private:
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:
    path_snapper(std::vector<point_t>& pts, std::vector<cell_t>& cells);
    path_snapper(simplicial::simplicial_complex&);
    ~path_snapper();
    path_snapper(path_snapper& other);
    path_snapper& operator=(path_snapper& other);

    std::vector<size_t> snap_path(std::vector<point_t> path);
    vector_t get_chain_vector(std::vector<point_t> path);
};
};
