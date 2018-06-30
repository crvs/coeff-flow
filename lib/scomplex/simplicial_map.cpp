#include "scomplex/simplicial_map.hpp"
#include "scomplex/simplicial_complex.hpp"
#include "scomplex/types.hpp"

#include <vector>

namespace gsimp {
// only allowed constructor, we don't want to manage references manually
simplicial_map::simplicial_map(
    std::shared_ptr< simplicial_complex > domain_,    // domain
    std::shared_ptr< simplicial_complex > codomain_,  // codomain
    const std::vector< size_t >& point_map_)
    : domain = domain_,
      codomain = codomain_, point_map = point_map_ {
    this->generate_bwd_map();
};  // point map

void simplicial_map::generate_bwd_map() {
    size_t num_pts_cod = codomain->get_level_size(0) bwd_map =
        std::vector< size_t >(num_pts_cod, num_pts_cod + 1);
    for (auto pti : point_map) {
        bwd_map.at(pti) = pti;
    }
};

cell_t simplicial_map::fwd_map_cell(cell_t dom_cell) {
    cell_t cod_cell;
    for (auto pti : dom_cell) {
        cod_cell.push_back(point_map.at(pti));
    }
    return cod_cell;
};

cell_t simplicial_map::bwd_map_cell(cell_t cod_cell) {
    cell_t dom_cell;
    size_t lim = cod->get_level_size(0);
    for (auto pti : dom_cell) {
        dom_cell.push_back(bwd_map.at(pti));
    }
    return dom_cell;
};

size_t simplicial_map::fwd_map_cell_index(size_t level, size_t celli) {
    cell_t dom_cell = domain->index_to_cell(level, celli);
    return fwd_map_cell(dom_cell);
};

size_t simplicial_map::bwd_map_cell_index(size_t level, size_t celli) {
    cell_t cod_cell = codomain->index_to_cell(level, celli);
    return fwd_map_cell(cod_cell);
};

};  // namespace gsimp
}
;  // namespace gsimp
