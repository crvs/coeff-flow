#include <vector>

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

#include "scomplex/coeff_flow.hpp"
#include "scomplex/simplicial_complex.hpp"
#include "scomplex/types.hpp"

using namespace gsimp;

namespace py = boost::python;

template < typename T >
inline std::vector< T > list_to_vector(const py::object& iterable) {
    return std::vector< T >(py::stl_input_iterator< T >(iterable),
                            py::stl_input_iterator< T >());
}

template < class T >
inline py::list vector_to_list(const std::vector< T >& v) {
    py::object get_iter = py::iterator< std::vector< T > >();
    py::object iter = get_iter(v);
    py::list l(iter);
    return l;
}

std::vector< cell_t > read_cells(const py::object& cells_obj) {
    std::vector< cell_t > cells_v;
    for (auto iter = py::stl_input_iterator< py::object >(cells_obj);  //
         iter != py::stl_input_iterator< py::object >();               //
         iter++) {                                                     //
        cells_v.push_back(list_to_vector< size_t >(*iter));
    }

    return cells_v;
}

class bounding_chain {
   public:
    std::shared_ptr< simplicial_complex > comp;

    bounding_chain(std::vector< cell_t > cells) {
        comp = std::make_shared< simplicial_complex >(cells);
    }

    size_t cell_to_index(cell_t cell) { return comp->cell_to_index(cell); }

    py::list get(std::vector< double > chain) {
        chain_v my_chain = {comp->dimension(), chain};
        chain_v bounding_chain = coeff_flow_embedded(*comp, my_chain);
        return vector_to_list< double >(bounding_chain.second);
    };
};

BOOST_PYTHON_MODULE(coeffflow) {
    using namespace boost::python;

    class_< bounding_chain >("bounding_chain",                 //
                             init< std::vector< cell_t > >())  //
        .def("cell_to_index", &bounding_chain::cell_to_index)  //
        .def("get", &bounding_chain::get);

    def("read_cells", read_cells);

    def("vector_to_list_d", vector_to_list<double>);
    def("vector_to_list_i", vector_to_list<int>);

}
