#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <vector>

int main() {
    struct Myoptions : Gudhi::Simplex_tree_options_full_featured {
        typedef int Vertex_handle;
    };
    Gudhi::Simplex_tree<Myoptions> simplexTree;
    auto d = {0, 1, 2};
    simplexTree.insert_simplex_and_subfaces(d);

    auto range = simplexTree.complex_simplex_range();
    int count = 0;

    for (auto s : range) {
        simplexTree.assign_key(s, count++);
        auto v_range = simplexTree.simplex_vertex_range(s);
        std::cout << simplexTree.filtration(s) << " , " << simplexTree.key(s)
                  << " : ";
        for (auto v : v_range) {
            std::cout << v << " ";
        }
        std::cout << ": ";
        auto sbd = simplexTree.boundary_simplex_range(s);
        for (auto sb : sbd) {
            auto v_range = simplexTree.simplex_vertex_range(sb);
            for (auto v : v_range) {
                std::cout << v << " ";
            }
            std::cout << " ; ";
        }
        std::cout << std::endl;
    }

    // testing the effects of reinsertion
    simplexTree.insert_simplex_and_subfaces(d);
    for (auto s : range) {
        auto v_range = simplexTree.simplex_vertex_range(s);
        std::cout << simplexTree.filtration(s) << " , " << simplexTree.key(s)
                  << " : ";
        for (auto v : v_range) {
            std::cout << v << " ";
        }
        std::cout << ": ";
        auto sbd = simplexTree.boundary_simplex_range(s);
        for (auto sb : sbd) {
            auto v_range = simplexTree.simplex_vertex_range(sb);
            for (auto v : v_range) {
                std::cout << v << " ";
            }
            std::cout << " ; ";
        }
        std::cout << std::endl;
    }

    auto p = simplexTree.find({0, 1, 2});
    for (auto v : simplexTree.simplex_vertex_range(p)) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::vector<int> vec;
    std::cout << "vector size " << vec.size() << std::endl;

    return 0;
}
