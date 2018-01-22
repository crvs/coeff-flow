#pragma once
#include "SimpComplex.hpp"
#include "chains.hpp"
#include "types.hpp"

#include <exception>
#include <tuple>
#include <vector>

#include <memory>
#include <fstream>

namespace Gsimp {

class NonBoundary : public std::exception {};
// NonBoundary is the exception thrown in case the given chain is not a boundary

/*
 * The BoundingChain object exists to automate the process of computing
 * bounding chains to boundaries on a given simplicial complex.
 *
 *
 * members:
 * s_comp:
 * - a shared_pointer the underlying simplicial complex (SimpComplex) object.
 * boundary_matrices:
 * - a vector containing at entry $i$ a (unique) pointer to the boundary matrix
 *   from dimension $i+1$ to dimension $i$.
 *
 * constructors:
 * BoundingChain(SimpComplex&)
 * BoundingChain(std::shared_ptr<SimpComplex>)
 * BoundingChain(std::vector<Point>&,std::vector<Cell>&)
 * - They instantiated the bounding chain object from the given data. In the
 *   first and second case, this is just computing the boundary matrices. In
 *   the third case the constructor also instantiates the underlying
 *   SimpComplex object.
 *
 * methods:
 * - populate_matrices(): create the boundary matrices of the simplicial complex
 * (needed before computing "getBoundingChain"
 * - getBoundingChain(chain) : takes a chain object and outputs a chain whose
 * boundary is the input chain (if it exists) and throws "NonBoundary"
 * exception otherwise.
 * */
class BoundingChain {
    std::shared_ptr< SimpComplex > s_comp;
    std::vector< std::unique_ptr< Matrix > > boundary_matrices;
    void populate_matrices();

   public:
    BoundingChain(SimpComplex&);
    BoundingChain(std::shared_ptr< SimpComplex >);
    BoundingChain(std::vector< Point >&, std::vector< Cell >&);
    ~BoundingChain();
    chain getBoundingChain(chain&);
};

};
