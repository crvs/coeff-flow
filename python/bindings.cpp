#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scomplex/path_snapper.hpp"
#include "scomplex/simplicial_complex.hpp"

#include "scomplex/coeff_flow.hpp"
#include "scomplex/types.hpp"

using namespace gsimp;

namespace py = pybind11;

class CoeffFlow {
    std::shared_ptr< simplicial_complex > base;
    std::shared_ptr< path_snapper > snapper;

   public:
    CoeffFlow(std::vector< cell_t >& tris) {
        base = std::make_shared< simplicial_complex >(tris);
        // snapper = std::make_shared< path_snapper >(base);
    }

    std::vector< double > bchain_from_paths(std::vector< point_t >& path1,
                                            std::vector< point_t >& path2) {
        chain_v ch1 = snapper->snap_path_to_v_chain(path1);
        chain_v ch2 = snapper->snap_path_to_v_chain(path2);
        for (size_t i = 0; i < std::get< 1 >(ch1).size(); i++) {
            std::get< 1 >(ch1).at(i) -= std::get< 1 >(ch2).at(i);
        }
        return std::get< 1 >(coeff_flow_embedded(*base, ch1));
    }

    std::vector< double > coeff_flow(std::vector< double >& chain) {
        return std::get< 1 >(coeff_flow_embedded(*base, {1, chain}));
    }

    std::vector< double > points_to_chain(std::vector< point_t >& pts) {
        return std::get< 1 >(snapper->snap_path_to_v_chain(pts));
    }
};

PYBIND11_MODULE(coeffflow, m) {
    py::class_< simplicial_complex >(m, "simplicial_complex")        //
        .def(py::init< std::vector< cell_t >& >())                   //
        .def("cell_to_index", &simplicial_complex::cell_to_index)    //
        .def("get_level_size", &simplicial_complex::get_level_size)  //
        .def("get_level", &simplicial_complex::get_level)            //
        .def("index_to_cell", &simplicial_complex::index_to_cell);

    py::class_< CoeffFlow >(m, "CoeffFlow")
        .def(py::init< std::vector< cell_t >& >())                //
        .def("coeff_flow", &CoeffFlow::coeff_flow)                //
        .def("bchain_from_paths", &CoeffFlow::bchain_from_paths)  //
        .def("points_to_chain", &CoeffFlow::points_to_chain);

    m.def("coeff_flow", coeff_flow);

    m.def("coeff_flow_embedded", coeff_flow_embedded);
};
