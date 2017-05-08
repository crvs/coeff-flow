#pragma once

#include <scomplex/types.hpp>
#include <scomplex/drawing_utils.hpp>
#include <scomplex/simplicial_complex.hpp>

namespace gsimp {

class DimensionException : public std::exception {};

void draw_complex(simplicial_complex& scomp, drawing& my_drawing,
                  double line_thickness, std::string fill_color = "#00ff00",
                  std::string line_color = "#000000") {
    if (scomp.dimension() > 2) DimensionException();

    auto points = scomp.get_points();
    auto edges = scomp.get_level(1);
    auto tris = scomp.get_level(2);

    for (auto cell : tris) {
        std::vector<std::vector<double>> cell_pt;
        for (auto i : cell) cell_pt.push_back(points[i]);
        my_drawing.draw_triangle(cell_pt, fill_color);
    }

    for (auto cell : scomp.get_level(1)) {
        std::vector<std::vector<double>> cell_pt;
        for (auto i : cell) cell_pt.push_back(points[i]);
        my_drawing.draw_edge(cell_pt[0], cell_pt[1], line_color,
                             line_thickness);
    }
}

void draw_path(path_snapper& p_snap, drawing& my_drawing,
               std::vector<std::vector<double>> path, double stroke,
               std::string path_color = "#ff0000") {
    auto snapped = p_snap.snap_path_to_points(path);
    auto it = snapped.begin();
    auto trail = it++;
    for (; it != snapped.end(); ++it, ++trail) {
        auto init = *trail;
        auto dest = *it;
        my_drawing.draw_edge(init, dest, path_color, stroke);
    }
}

void draw_2_chain(simplicial_complex& s_comp, drawing& my_drawing,
                  chain_t chain, std::string chain_color = "#ff0000") {
    if (chain_dim(chain) != 2) DimensionException();

    auto points = s_comp.get_points();
    auto level_2 = s_comp.get_level(2);
    Eigen::SparseVector<double>::InnerIterator it(chain_rep(chain));
    for (; it; ++it) {
        cell_t cell = level_2[it.index()];
        std::vector<point_t> cell_pt;
        for (size_t i : cell) cell_pt.push_back(points[i]);
        // TODO: get rid of these thresholds
        if (it.value() > 0.01 || it.value() < -0.01)
        my_drawing.draw_triangle(cell_pt, chain_color);
    }
}
void draw_1_chain(simplicial_complex& s_comp, drawing& my_drawing,
                  chain_t chain, double stroke,
                  std::string chain_color = "#ff0000") {
    if (chain_dim(chain) != 1) DimensionException();

    auto points = s_comp.get_points();
    auto level_1 = s_comp.get_level(1);
    Eigen::SparseVector<double>::InnerIterator it(chain_rep(chain));
    for (; it; ++it) {
        cell_t edge = level_1[it.index()];
        my_drawing.draw_edge(points[edge[0]], points[edge[1]], chain_color,
                             stroke);
    }
}

void draw_path(simplicial_complex& scomp, drawing& my_drawing,
               std::vector<std::vector<double>> path, double stroke,
               std::string path_color = "#ff0000") {
    path_snapper p_snap(scomp);
    auto snapped = p_snap.snap_path_to_points(path);
    auto it = snapped.begin();
    auto trail = it++;
    for (; it != snapped.end(); ++it, ++trail) {
        auto init = *trail;
        auto dest = *it;
        my_drawing.draw_edge(init, dest, path_color, stroke);
    }
}

drawing draw_complex(simplicial_complex& scomp, std::string filename,
                     double stroke, int width = 500, int height = 500,
                     std::string fill_color = "#00ff00",
                     std::string line_color = "#000000") {
    drawing my_drawing(filename, width, height);
    draw_complex(scomp, my_drawing, stroke, fill_color, line_color);
    return my_drawing;
}
};
