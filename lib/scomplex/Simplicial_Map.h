#ifndef SIMPLICIAL_MAP
#define SIMPLICIAL_MAP

#include <vector>

#include <scomplex/Simplicial_Complex.h>
namespace simplicial {

class SimplicialMap {
    SimplicialComplex* domain;
    SimplicialComplex* codomain;

    // type of index transformations
    typedef std::vector<int> ind_trans;

    // transformation of points
    ind_trans point_trans;

    // the map between the simplices
    std::vector<ind_trans> simplicial_map;
    // simplicial_map[i] -- is an ind_trans between the indices of the domain
    // and the codomain

    chain_t image(chain_t chain) { return chain; };

    chain_t preimage(chain_t chain) { return chain; };

};      // class SimplicialMap
};      // namespace simplicial
#endif  // SIMPLICIAL_MAP
