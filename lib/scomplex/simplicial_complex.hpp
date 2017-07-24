#pragma once

#include <scomplex/types.hpp>
#include <scomplex/chains.hpp>
#include <memory>

#include <gudhi/Simplex_tree.h>

namespace gsimp {


class No_Boundary {};

class simplicial_complex {
    // implementation details
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:

    void calculate_hasse();

    // constructor (no default)
    simplicial_complex(const std::vector<cell_t>&);
    simplicial_complex(const std::vector<point_t>&,const std::vector<cell_t>&);
    simplicial_complex(const simplicial_complex&);
    simplicial_complex& operator=(const simplicial_complex&);
    // destructor
    ~simplicial_complex();
    // basic info
    std::vector<point_t> get_points();
    point_t get_point(size_t);
    int dimension();
    // level-wise info
    chain new_sparse_chain(int d);
    chain new_dense_chain(int d);
    chain new_chain(int d);
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
    // boundary matrices
    matrix_t get_boundary_matrix(int);
    // cells and indices back and forth
    cell_t index_to_cell(int, size_t);
    size_t cell_to_index(cell_t);
    // TODO
    // volume_chain(int d) {
    // /* returns a cahin containing the volumes of all d--cells of the complex
    //  * computes it once stores it in pimpl and retrieves it whenever it becomes necessary
    //  * it returns a dense chain
    //  */
    // }
};  // class simplicial_complex
};  // namespace gsimp
