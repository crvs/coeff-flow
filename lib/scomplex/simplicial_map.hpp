#pragma once

#include "scomplex/simplicial_complex.hpp"
#include "scomplex/types.hpp"

#include <vector>

namespace gsimp {
class simplicial_map {
   private:
    std::shared_ptr< simplicial_complex > domain;
    std::shared_ptr< simplicial_complex > codomain;

    const std::vector< size_t > fwd_map;
    const std::vector< size_t > bwd_map;

    void generate_bwd_map();

   public:
    // only allowed constructor, we don't want to manage references manually
    simplicial_map(std::shared_ptr< simplicial_complex >,  // domain
                   std::shared_ptr< simplicial_complex >,  // codomain
                   const std::vector< size_t >&);           // point map

    cell_t fwd_map_cell(cell_t);
    cell_t bwd_map_cell(cell_t);

    size_t fwd_map_cell_index(size_t, size_t);
    size_t bwd_map_cell_index(size_t, size_t);

};  // class simplicial_map
};  // namespace gsimp
