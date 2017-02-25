#pragma once

#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <cmath>
#include <exception>
#include <tuple>
#include <vector>

#include <memory>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <fstream>

namespace gsimp {

class non_zero_chain : public std::exception {};
class cant_draw_dimension : public std::exception {};

class bounding_chain {
    std::shared_ptr<simplicial_complex> s_comp;
    std::vector<std::unique_ptr<matrix_t>> boundary_matrices;
    void populate_matrices();

   public:
    bounding_chain(simplicial_complex& sc);
    bounding_chain(std::shared_ptr<simplicial_complex> sc);
    bounding_chain(std::vector<point_t>& points, std::vector<cell_t>& tris);
    ~bounding_chain();
    chain_t get_bounding_chain(chain_t&);
    // void draw_chain(std::string, chain_t);
};

//---------------------------------------
// implementation

bounding_chain::~bounding_chain() {}

bounding_chain::bounding_chain(std::shared_ptr<simplicial_complex> sc) {
    s_comp = sc;
    populate_matrices();
}

bounding_chain::bounding_chain(simplicial_complex& sc) {
    s_comp.reset(new simplicial_complex(sc));
    populate_matrices();
}

bounding_chain::bounding_chain(std::vector<point_t>& points,
                               std::vector<cell_t>& tris) {
    s_comp.reset(new simplicial_complex(points, tris));
    populate_matrices();
}

void bounding_chain::populate_matrices() {
    for (int d = 0; d < s_comp->dimension(); ++d) {
        auto level_matrix = s_comp->get_boundary_matrix(d);
        boundary_matrices.push_back(
            std::unique_ptr<matrix_t>(new matrix_t(level_matrix)));
    }
}
//
// round vectors for comparison
vector_t round_vec(vector_t vec) {
    vector_t rounded_vec(vec.rows());
    for (int i = 0; i < vec.rows(); ++i) {
        int coef = round(vec.coeffRef(i));
        if (coef != 0) rounded_vec.coeffRef(i) = coef;
    }
    return rounded_vec;
}

bool equals(vector_t vec1, vector_t vec2) {
    for (int i = 0; i < vec1.rows(); ++i) {
        if (vec1.coeffRef(i) != vec2.coeffRef(i)) {
            std::cout << i << ": " << vec1.coeffRef(i) << " "
                      << vec2.coeffRef(i) << '\n';
            return false;
        }
    }
    return true;
}

// // TODO: drowing facilities
// //  - draw chain
// //  - draw 1-chain with arrows
// //  - put it into its own header file
// void bounding_chain::draw_chain(std::string filename, chain_t chain) {
//     typedef boost::geometry::model::d2::point_xy<double> point_type;
//     typedef boost::geometry::model::ring<point_type> polygon_type;
//
//     std::vector<polygon_type> polys;
//
//     for (auto tri : s_comp->get_level(2)) {
//         std::vector<point_type> polygon_vec;
//         for (auto k : tri) {
//             point_t pt = s_comp->get_points().at(k);
//             point_type pt_g(pt[0], pt[1]);
//             polygon_vec.push_back(pt_g);
//         }
//         polygon_type poly(polygon_vec.begin(), polygon_vec.end());
//         polys.push_back(poly);
//     }
//
//     // now we want to draw this
//     std::ofstream svg(filename.c_str());
//
//     boost::geometry::svg_mapper<point_type> mapper(svg, 150, 150);
//
//     int chain_d = std::get<0>(chain);
//     vector_t chain_v = std::get<1>(chain);
//     for (int i = 0; i < chain_v.rows(); ++i) {
//         polygon_type poly = polys.at(i);
//         mapper.add(poly);
//         if (chain_d == 2) {
//             if (chain_v.coeffRef(i) == 0)
//                 mapper.map(
//                     poly,
//                     "fill-opacity:0.5;fill:rgb(153,204,0);stroke:rgb(153,"
//                     "204,0);stroke-width:2");
//             else
//                 mapper.map(
//                     poly,
//                     "fill-opacity:0.5;fill:rgb(230,104,58);stroke:rgb(230,"
//                     "104,58);stroke-width:2");
//         } else {
//             mapper.map(poly,
//                        "fill-opacity:0.5;fill:rgb(153,204,0);stroke:rgb(153,"
//                        "204,0);stroke-width:2");
//         }
//
//         if (chain_d == 1) {
//             auto lines = s_comp->get_level(1);
//             for (int i = 0; i < lines.size(); ++i) {
//                 if (chain_v.coeffRef(i) != 0) {
//                     std::vector<point_type> line_vec;
//                     for (auto k : lines.at(i)) {
//                         point_t pt = s_comp->get_points().at(k);
//                         point_type pt_g(pt[0], pt[1]);
//                         line_vec.push_back(pt_g);
//                     }
//                     polygon_type line_g(line_vec.begin(), line_vec.end());
//                     mapper.map(
//                         line_g,
//                         "fill-opacity:0.5;fill:rgb(250,0,0);stroke:rgb(250,"
//                         "0,0);stroke-width:2");
//                 }
//             }
//         }
//     }
// }

chain_t bounding_chain::get_bounding_chain(chain_t& chain) {
    int chain_d;
    vector_t chain_v;
    std::tie<int, vector_t>(chain_d, chain_v) = chain;

    if (chain_d >= s_comp->dimension()) throw non_zero_chain();

    Eigen::LeastSquaresConjugateGradient<matrix_t> lscg;
    lscg.compute(*(boundary_matrices.at(chain_d)));

    vector_t bound_chain(round_vec(lscg.solve(chain_v)));
    vector_t result(round_vec(*(boundary_matrices.at(chain_d)) * bound_chain));

    if (equals(result, chain_v))
        return chain_t(chain_d - 1, bound_chain);
    else {
        // std::cout << "where to write the svg file?";
        // std::string filename;
        // std::cin >> filename;
        // draw_chain(filename, chain_t(1, bound_chain));
        // std::cout << '\n';
        // std::cout << result << '\n';
        // std::cout << chain_v << '\n';
        throw non_zero_chain();
    }
}
};
