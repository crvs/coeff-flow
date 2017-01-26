#include <scomplex/Simplicial_Complex.h>
#include <scomplex/utils.h>
#include <scomplex/trace.h>

#include "boost/graph/graph_traits.hpp"

#include <iostream>
#include <vector>
#include <string>

typedef std::vector<double> point_t;
typedef std::list<int> simp_t;
typedef std::pair<int, std::vector<int>> chain_t;

int f(point_t a) { return a.at(0) < .5 ? 0 : 1; }

int main() {
    std::list<point_t> point_list;
    std::list<simp_t> simp_list;
    // ---

    point_list.push_back(point_t({0, 0}));
    point_list.push_back(point_t({0, 1}));
    point_list.push_back(point_t({1, 0}));
    point_list.push_back(point_t({1, 1}));

    simp_list.push_back(simp_t({0, 1, 2}));
    simp_list.push_back(simp_t({1, 2, 3}));

    simplicial::SimplicialComplex<point_t> sc(point_list, simp_list);
    for (auto pt : sc.points) {
        std::cout << pt.at(0) << " " << pt.at(1) << std::endl;
    }

    auto range = sc.simplices.complex_simplex_range();
    for (auto s : range) {
        std::cout << sc.simplices.dimension(s) << " , " << sc.simplices.key(s)
                  << " : ";
        for (auto v : sc.simplices.simplex_vertex_range(s)) {
            std::cout << v << " ";
        }
        // iterate over boundary
        std::cout << std::endl
                  << "    boundary: ";
        for (auto sb : sc.simplices.boundary_simplex_range(s)) {
            auto v_range = sc.simplices.simplex_vertex_range(sb);
            for (auto v : v_range) {
                std::cout << v << " ";
            }
            std::cout << " ; ";
        }
        std::cout << std::endl
                  << "  coboundary: ";
        for (auto sb : sc.simplices.cofaces_simplex_range(s, 1)) {
            auto v_range = sc.simplices.simplex_vertex_range(sb);
            for (auto v : v_range) {
                std::cout << v << " ";
            }
            std::cout << " ; ";
        }
        std::cout << std::endl;
    }

    std::cout << "now quotient" << std::endl;
    auto sq = sc.quotient(f);

    /* the points are working, what's up with the simplices??
    for (auto p : sq.points) {
        std::cout << "point: " << p.at(0) << " " << p.at(1) << std::endl;
    }
    */

    range = sq.simplices.complex_simplex_range();
    std::cout << "got this far" << std::endl;
    for (auto s : range) {
        std::cout << sq.simplices.dimension(s) << " , " << sq.simplices.key(s)
                  << " : ";
        for (auto v : sq.simplices.simplex_vertex_range(s)) {
            std::cout << v << " ";
        }
        // iterate over boundary
        std::cout << std::endl
                  << "    boundary: ";
        for (auto sb : sq.simplices.boundary_simplex_range(s)) {
            auto v_range = sq.simplices.simplex_vertex_range(sb);
            for (auto v : v_range) {
                std::cout << v << " ";
            }
            std::cout << " ; ";
        }
        std::cout << std::endl
                  << "  coboundary: ";
        for (auto sb : sq.simplices.cofaces_simplex_range(s, 1)) {
            auto v_range = sq.simplices.simplex_vertex_range(sb);
            for (auto v : v_range) {
                std::cout << v << " ";
            }
            std::cout << " ; ";
        }
        std::cout << std::endl;
    }
    std::cout << "got this far" << std::endl;

    // checking the matrices
    auto g = sc.get_boundary(1);
    std::cout << g << std::endl;
    g = sc.get_boundary(0);
    std::cout << g << std::endl;

    g = sq.get_boundary(1);
    std::cout << g << std::endl;
    g = sq.get_boundary(0);
    std::cout << g << std::endl;

    std::cout << "printing level: " << std::endl;
    for (auto s : sq.get_level(1)) {
        std::cout << s.size() << ": ";
        std::cout << s.at(0) << " ";
        std::cout << s.at(1) << " ";
        std::cout << std::endl;
    }

    auto G = calculate_one_skelleton_graph(sc);
    auto zero = Vertex(0, G);
    auto one = Vertex(1, G);
    auto shortest_path(G, zero, one)

        return 0;
}
