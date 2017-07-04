#pragma once

#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>

namespace gsimp {

class path_snapper {
   private:
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:
    // constructors and such
    path_snapper(std::vector<point_t>&, std::vector<cell_t>&);
    path_snapper(std::shared_ptr<simplicial_complex>);
    path_snapper(simplicial_complex&);
    ~path_snapper();
    path_snapper(path_snapper&);
    path_snapper& operator=(const path_snapper&);
    // the things we want to do
    std::vector<point_t> snap_path_to_points(std::vector<point_t>);
    std::vector<size_t> snap_path_to_indices(std::vector<point_t>);
    chain_t snap_path_to_chain(std::vector<point_t>);
    chain_v snap_path_to_v_chain(std::vector<point_t>);
    // interconversion
    std::vector<point_t> index_sequence_to_point(std::vector<size_t>);
    std::vector<size_t> point_sequence_to_index(std::vector<point_t>);
    chain_t index_sequence_to_chain(std::vector<size_t>);
    chain_v index_sequence_to_v_chain(std::vector<size_t>);
    chain_t point_sequence_to_chain(std::vector<point_t>);
    std::shared_ptr<simplicial_complex> get_underlying_complex();
};
};
