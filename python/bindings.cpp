#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scomplex/coeff_flow.hpp"
#include "scomplex/simplicial_complex.hpp"
#include "scomplex/types.hpp"

using namespace gsimp;

namespace py = pybind11;

PYBIND11_MODULE(coeffflow, m) {
    py::class_< simplicial_complex >(m, "simplicial_complex")  //
        .def(py::init< std::vector< cell_t >& >())             //
        .def("cell_to_index", &simplicial_complex::cell_to_index);

    m.def("coeff_flow", coeff_flow);

    m.def("coeff_flow_embedded",coeff_flow_embedded);
};
