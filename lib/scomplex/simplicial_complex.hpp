#pragma once

#include <scomplex/types.hpp>
#include <memory>

#include <gudhi/Simplex_tree.h>

namespace gsimp {


class NoChain {};

class simplicial_complex {
    // implementation details
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:

    void calculate_hasse();

    // constructor (no default)
    simplicial_complex(std::vector<point_t>&, std::vector<cell_t>&);
    simplicial_complex(const simplicial_complex&);
    simplicial_complex& operator=(const simplicial_complex&);
    // destructor
    ~simplicial_complex();
    // basic info
    std::vector<point_t> get_points();
    point_t get_point(size_t);
    int dimension();
    bool is_quotient();
    // level-wise info
    chain_t new_chain(int d);
    chain_v new_v_chain(int d);
    // create chains
    std::vector<cell_t> get_level(int);
    int get_level_size(int);
    // inclusion index
    int boundary_inclusion_index(cell_t, cell_t);
    int boundary_inclusion_index(int, size_t, int, size_t);
    // cell boundaries
    std::vector<cell_t> cell_boundary(cell_t);
    std::vector<size_t> cell_boundary_index(int, size_t);
    std::vector<std::pair<int, cell_t>> get_bdry_and_ind(cell_t);
    std::vector<std::pair<int, size_t>> get_bdry_and_ind_index(int, size_t);
    // treating cofaces
    std::vector<cell_t> get_cofaces(cell_t);
    std::vector<size_t> get_cofaces_index(int, size_t);
    std::vector<std::pair<int, cell_t>> get_cof_and_ind(cell_t);
    std::vector<std::pair<int, size_t>> get_cof_and_ind_index(int, size_t);
    // cell area calculations
    double cell_area(cell_t);
    double cell_index_area(int,size_t);
    double chain_area(chain_t);
    // boundary matrices
    matrix_t get_boundary_matrix(int);
    // cells and indices back and forth
    cell_t index_to_cell(int, size_t);
    size_t cell_to_index(cell_t);
    //
    simplicial_complex quotient(int char_fun(point_t));
};  // class simplicial_complex
};  // namespace gsimp
