#include <scomplex/Simplicial_Complex.h>

#include <iostream>
#include <vector>

typedef std::vector<float> point_t;
typedef std::list<int> simp_t;
typedef std::pair<int, std::vector<int>> chain_t;

int f(point_t a) {
    int b;
    b = a.at(0) < .5 ? 0 : 1;
    return b;
}

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

    auto sc = simplicial::SimplicialComplex<point_t>(point_list, simp_list);

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

    auto sq = sc.quotient(f);
    std::cout << std::endl;
    std::cout << "now quotient" << std::endl;
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

    // checking the matrices
    auto g = sc.get_boundary(1);
    std::cout << g << std::endl;
    g = sc.get_boundary(0);
    std::cout << g << std::endl;

    g = sq.get_boundary(1);
    std::cout << g << std::endl;
    g = sq.get_boundary(0);
    std::cout << g << std::endl;

    return 0;
}
