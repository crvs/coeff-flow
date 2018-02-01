#pragma once

#include <boost/function.hpp>

#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>

/*
 *
 */

namespace gsimp {
class quotient{

    struct impl;
    std::unique_ptr<impl> p_impl;

    public:
        quotient(const simplicial_complex&, boost::function<bool(point_t)>);
        quotient(const simplicial_complex&, boost::function<bool(point_t)>, point_t);
        quotient(std::shared_ptr<simplicial_complex>, boost::function<bool(point_t)>);
        quotient(std::shared_ptr<simplicial_complex>, boost::function<bool(point_t)>, point_t);
        quotient(const quotient&);
        quotient& operator=(const quotient&);

        ~quotient();

        std::shared_ptr<simplicial_complex> base_complex();
        std::shared_ptr<simplicial_complex> quotient_complex();

        // push chains forward through quotient
        chain quotient_chain(chain);

        // pull chains backward through quotient
        chain unquotient_chain(chain);

}; // quotient
}; // gsimp
