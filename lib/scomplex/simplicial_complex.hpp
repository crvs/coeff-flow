#pragma once

#include <scomplex/types.hpp>
#include <memory>

namespace gsimp {

using namespace gsimp;

class NoChain {};

class simplicial_complex {
    // implementation details
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:
    // constructor (no default)
    simplicial_complex(std::vector<point_t>&, std::vector<cell_t>&);
    simplicial_complex(const simplicial_complex&);
    simplicial_complex& operator=(const simplicial_complex&);
    // destructor
    ~simplicial_complex();
    //
    std::vector<point_t> get_points();
    std::vector<cell_t> get_level(int);
    int dimension();
    int get_level_size(int);
    matrix_t get_boundary_matrix(int);
    bool is_quotient();
    size_t get_index_of_simplex(cell_t);
    //
    simplicial_complex quotient(int char_fun(point_t));
};  // class simplicial_complex
};  // namespace gsimp
