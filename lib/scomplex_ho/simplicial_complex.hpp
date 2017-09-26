#pragma once

#include <scomplex_ho/exterior_prod.hpp>
#include <scomplex_ho/chains.hpp>
#include <scomplex_ho/types.hpp>

#include <iterator>  // for debuging purposes

#include <gudhi/Simplex_tree.h>
#include <gudhi/Hasse_complex.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <tuple>
#include <memory>  // smart pointers
#include <functional>
#include <cmath>

#include <iostream>

namespace gsimp {

class No_Boundary {};

class simplicial_complex {
    // implementation details
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:

    void calculate_hasse();

    // constructor (no default)
    simplicial_complex(const std::vector<cell_t>&);
    simplicial_complex(const std::vector<point_t>&,const std::vector<cell_t>&);
    simplicial_complex(const simplicial_complex&);
    simplicial_complex& operator=(const simplicial_complex&);
    // destructor
    ~simplicial_complex();
    // basic info
    std::vector<point_t> get_points();
    point_t get_point(size_t);
    int dimension();
    // level-wise info
    chain new_sparse_chain(int d);
    chain new_dense_chain(int d);
    chain new_chain(int d);
    // create chains
    std::vector<cell_t> get_level(int);
    int get_level_size(int);
    // inclusion index
    int boundary_inclusion_index(cell_t, cell_t);
    int boundary_inclusion_index(int, size_t, int, size_t);
    // cell boundaries
    std::vector<cell_t> cell_boundary(cell_t);
    std::vector<size_t> cell_boundary_index(int, size_t);
    std::vector<std::pair<int, cell_t>> get_bdry_and_ind(cell_t);
    std::vector<std::pair<int, size_t>> get_bdry_and_ind_index(int, size_t);
    // treating cofaces
    std::vector<cell_t> get_cofaces(cell_t);
    std::vector<size_t> get_cofaces_index(int, size_t);
    std::vector<std::pair<int, cell_t>> get_cof_and_ind(cell_t);
    std::vector<std::pair<int, size_t>> get_cof_and_ind_index(int, size_t);
    // boundary matrices
    matrix_t get_boundary_matrix(int);
    // cells and indices back and forth
    cell_t index_to_cell(int, size_t);
    size_t cell_to_index(cell_t);
    chain volume_chain(int);
    // /* returns a chain containing the volumes of all d--cells of the complex
    //  * computes it once stores it in pimpl and retrieves it whenever it becomes necessary
    //  * it returns a dense chain
    //  */
    // }
};  // class simplicial_complex

// hasse diagram structures
struct hasse_node {
    std::pair<int, size_t> handle;
    std::vector<std::shared_ptr<hasse_node>> cofaces;

    hasse_node(std::pair<int, size_t> _handle) { handle = _handle; }

    void add_coface(std::shared_ptr<hasse_node>& node) {
        cofaces.push_back(node);
    }
};

struct hasse_diag {
    std::vector<std::vector<std::shared_ptr<hasse_node>>> cells;
    std::vector<std::vector<bool>> cells_c;

    hasse_diag(){};

    hasse_diag(simplicial_complex& s_comp) {
        // vector of dimension d
        for (int d = 0; d <= s_comp.dimension(); d++) {
            size_t lv_s = s_comp.get_level_size(d);
            std::vector<std::shared_ptr<hasse_node>> lv(lv_s);
            cells.push_back(lv);

            std::vector<bool> c(lv_s, false);
            cells_c.push_back(c);
        }

        for (int d = s_comp.dimension(); d > 0; d--) {
            for (size_t face_i = 0; face_i < s_comp.get_level_size(d);
                 face_i++) {
                std::shared_ptr<hasse_node> face_ptr =
                    get_face(d, face_i, true);

                std::vector<std::pair<int, size_t>> faces =
                    s_comp.get_bdry_and_ind_index(d, face_i);
                for (auto b_face_i : faces) {
                    size_t b_face = std::get<1>(b_face_i);
                    std::shared_ptr<hasse_node> b_face_ptr =
                        get_face(d - 1, b_face, true);
                    b_face_ptr->add_coface(face_ptr);
                }
            }
        }
        // no further need, free up the space
        cells_c.clear();
    }

    std::shared_ptr<hasse_node> get_face(int d, size_t face_i,
                                         bool building = false) {
        std::shared_ptr<hasse_node> face;
        if (building && cells_c[d][face_i])
            face = cells[d][face_i];
        else {
            if (building) cells_c[d][face_i] = true;

            cells[d][face_i] = std::shared_ptr<hasse_node>(
                new hasse_node(std::pair<int, size_t>(d, face_i)));
            face = cells[d][face_i];
        }
        return face;
    }

    std::vector<size_t> get_coface_i(int d, size_t face_i) {
        auto node = cells[d][face_i];
        std::vector<size_t> cofaces;
        for (auto v : node->cofaces) cofaces.push_back(std::get<1>(v->handle));
        return cofaces;
    }
};

struct simplicial_complex::impl {
    // auxiliary types
    struct SimpleOptions : Gudhi::Simplex_tree_options_full_featured {
        typedef size_t Vertex_handle;
    };

    typedef Gudhi::Simplex_tree<SimpleOptions> simp_tree;
    typedef Gudhi::Simplex_tree<SimpleOptions>::Simplex_handle simp_handle;
    typedef std::vector<simp_handle*> level_t;
    typedef std::vector<std::unique_ptr<level_t>> levels_t;

    // member variables
    std::vector<point_t> points;
    simp_tree simplices;
    std::vector<matrix_t> boundary_matrices;
    levels_t levels;

    bool has_hasse;
    hasse_diag incidence;

    impl(const std::vector<point_t>& arg_points,
         const std::vector<cell_t>& arg_tris)
        : points(arg_points) {
        // create the simplex tree
        has_hasse=false;
        for (auto tri : arg_tris) {
            // removed deduping to try to make this a bit faster
            // ... it did cut time down about 10%, so ...
            simplices.insert_simplex_and_subfaces(tri);
            int d = tri.size() - 1;
            if (simplices.dimension() < d) simplices.set_dimension(d);
        }

        // assign a key to each simplex in each level
        std::vector<size_t> count(simplices.dimension() + 1, 0);
        for (int i = 0; i < simplices.dimension() + 1; ++i) {
            level_t* level = new level_t();
            levels.push_back(std::unique_ptr<level_t>(level));
        }

        auto simplex_range = simplices.complex_simplex_range();
        for (auto s : simplex_range) {
            int d = simplices.dimension(s);
            simplices.assign_key(s, count[d]++);
            levels[d]->push_back(new simp_handle(s));
        }
    }

    ~impl() {
        // get rid of the levels so they won't dangle
        for (int i = levels.size(); i <= 0; --i) levels[i].reset();
    };

    size_t get_level_size(int level) { return level<levels.size() ? levels[level]->size() : 0; }

    // calculate the index of s_1 in the boundary of s_2
    int boundary_index(simp_handle s_1, simp_handle s_2) {
        auto it_1 = simplices.simplex_vertex_range(s_1).begin();
        auto end_1 = simplices.simplex_vertex_range(s_1).end();
        //
        auto it_2 = simplices.simplex_vertex_range(s_2).begin();
        //
        int orient = 0;
        int unmatch = 0;
        //
        for (int p = 0; p <= simplices.dimension(s_2); p++) {
            if (*it_1 != *it_2) {
                orient = p++;
                unmatch++;
            } else {
                if (it_1 != end_1) it_1++;
                it_2++;
            }
        }
        return pow(-1, orient);
    }

    std::vector<size_t> dedupe_vec(std::vector<size_t>& vec) {
        std::set<size_t> no_reps;
        std::vector<size_t> no_reps_list;
        for (size_t el : vec) {
            no_reps.insert(el);
        }
        for (size_t el : no_reps) {
            no_reps_list.push_back(el);
        }
        return no_reps_list;
    }

    void calculate_matrices() {
        boundary_matrices = std::vector<matrix_t>();
        for (int k = 0; k < simplices.dimension(); k++) {
            boundary_matrices.push_back(
                matrix_t(get_level_size(k), get_level_size(k + 1)));
        }
        for (auto s : simplices.complex_simplex_range()) {
            int j = simplices.key(s);
            for (auto bs : simplices.boundary_simplex_range(s)) {
                int i = simplices.key(bs);
                int k = simplices.dimension(bs);
                boundary_matrices[k].coeffRef(i, j) = boundary_index(bs, s);
            }
        }
    }

    std::vector<cell_t> get_level(int level) {
        std::vector<cell_t> level_cells;
        for (auto simp : *levels[level]) {
            cell_t v_simp;
            for (auto v : simplices.simplex_vertex_range(*simp)) {
                v_simp.push_back(v);
            }
            level_cells.push_back(v_simp);
        }
        return level_cells;
    }

    simp_handle index_to_handle(int d, size_t tau) {
        return *(*(levels[d]))[tau];
    }

    size_t handle_to_index(simp_handle tau) { return simplices.key(tau); }

    simp_handle cell_to_handle(cell_t tau) {
        auto sh = simplices.find(tau);
        return sh;
    }

    cell_t handle_to_cell(simp_handle tau) {
        cell_t cell;
        for (auto v : simplices.simplex_vertex_range(tau)) cell.push_back(v);
        return cell;
    }

};  // struct impl

std::vector<std::pair<int, cell_t>> simplicial_complex::get_bdry_and_ind(
    cell_t cell) {
    impl::simp_handle simp = p_impl->cell_to_handle(cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<std::pair<int, cell_t>> boundary_and_indices;
    for (auto face : c_boundary)
        boundary_and_indices.push_back(              //
            std::make_pair<int, cell_t>(             //
                p_impl->boundary_index(face, simp),  //
                p_impl->handle_to_cell(face)));      //
    return boundary_and_indices;
};

std::vector<std::pair<int, size_t>> simplicial_complex::get_bdry_and_ind_index(
    int d, size_t cell) {
    impl::simp_handle simp = p_impl->index_to_handle(d, cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<std::pair<int, size_t>> boundary_and_indices;
    for (auto face : c_boundary)
        boundary_and_indices.push_back(              //
            std::make_pair<int, size_t>(             //
                p_impl->boundary_index(face, simp),  //
                p_impl->handle_to_index(face)));     //
    return boundary_and_indices;
};

std::vector<size_t> simplicial_complex::cell_boundary_index(int d,
                                                            size_t cell) {
    impl::simp_handle simp = p_impl->index_to_handle(d, cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<size_t> s_boundary;
    for (auto c : c_boundary) s_boundary.push_back(p_impl->handle_to_index(c));
    return s_boundary;
}

std::vector<cell_t> simplicial_complex::cell_boundary(cell_t cell) {
    impl::simp_handle simp = p_impl->cell_to_handle(cell);
    auto c_boundary = p_impl->simplices.boundary_simplex_range(simp);
    std::vector<cell_t> s_boundary;
    for (auto c : c_boundary) s_boundary.push_back(p_impl->handle_to_cell(c));
    return s_boundary;
}

int simplicial_complex::boundary_inclusion_index(cell_t c1, cell_t c2) {
    return p_impl->boundary_index(p_impl->cell_to_handle(c1),   //
                                  p_impl->cell_to_handle(c2));  //
};
int simplicial_complex::boundary_inclusion_index(int d1, size_t s1,    //
                                                 int d2, size_t s2) {  //
    cell_t c1 = index_to_cell(d1, s1);
    cell_t c2 = index_to_cell(d2, s2);
    return p_impl->boundary_index(p_impl->cell_to_handle(c1),
                                  p_impl->cell_to_handle(c2));
};

std::vector<std::pair<int, cell_t>> simplicial_complex::get_cof_and_ind(
    cell_t cell) {
    std::vector<std::pair<int, cell_t>> c_cofaces;
    for (auto face : get_cofaces(cell))
        c_cofaces.push_back(                                   //
            std::make_pair(                                    //
                boundary_inclusion_index(cell, face), face));  //
    return c_cofaces;
}

std::vector<std::pair<int, size_t>> simplicial_complex::get_cof_and_ind_index(
    int d, size_t c) {
    std::vector<std::pair<int, size_t>> c_cofaces;
    for (auto face : get_cofaces_index(d, c))
        c_cofaces.push_back(                                          //
            std::make_pair(                                           //
                boundary_inclusion_index(d, c, d + 1, face), face));  //
    return c_cofaces;
}

int simplicial_complex::get_level_size(int level) {
    return p_impl->get_level_size(level);
}

simplicial_complex::simplicial_complex(const std::vector<point_t>& arg_points,
                                       const std::vector<cell_t>& arg_tris) {
    p_impl = std::shared_ptr<impl>(new impl(arg_points, arg_tris));
}

simplicial_complex::simplicial_complex(const std::vector<cell_t>& arg_tris) {
    std::vector<point_t> arg_points;
    p_impl = std::shared_ptr<impl>(new impl(arg_points, arg_tris));
}

simplicial_complex::simplicial_complex(const simplicial_complex& other) {
    p_impl = other.p_impl;
}

simplicial_complex& simplicial_complex::operator=(
    const simplicial_complex& other) {
    p_impl = other.p_impl;
    return *this;
}

simplicial_complex::~simplicial_complex() {}

std::vector<point_t> simplicial_complex::get_points() { return p_impl->points; }

point_t simplicial_complex::get_point(size_t index) {
    return p_impl->points[index];
}

std::vector<cell_t> simplicial_complex::get_level(int level) {
    return p_impl->get_level(level);
}

matrix_t simplicial_complex::get_boundary_matrix(int d) {
    // uninstantiated boundary matrices
    if (p_impl->boundary_matrices.size() == 0) p_impl->calculate_matrices();

    // now they have to be instantiated, get them
    if (0 <= d && d < p_impl->boundary_matrices.size())
        return p_impl->boundary_matrices[d];
    else
        throw No_Boundary();
}

int simplicial_complex::dimension() { return p_impl->simplices.dimension(); }

cell_t simplicial_complex::index_to_cell(int d, size_t ind) {
    auto sh = (*(p_impl->levels[d]))[ind];
    return p_impl->handle_to_cell(*sh);
}

size_t simplicial_complex::cell_to_index(cell_t simp) {
    auto sh = p_impl->simplices.find(simp);
    return p_impl->simplices.key(sh);
}

std::vector<size_t> simplicial_complex::get_cofaces_index(int d, size_t face) {
    // codimension 1 faces
    if (!p_impl->has_hasse) {
        calculate_hasse();
    }
    auto s_cofaces = p_impl->incidence.get_coface_i(d, face);
    return s_cofaces;
}

std::vector<cell_t> simplicial_complex::get_cofaces(cell_t face) {
    if (!p_impl->has_hasse) {
        calculate_hasse();
    }
    std::vector<cell_t> s_cofaces;
    auto face_i = cell_to_index(face);
    int d = face.size() - 1;
    // codimension 1 faces
    auto coface_i_v = p_impl->incidence.get_coface_i(d, face_i);
    for (auto v : coface_i_v) {
        auto face_h = (*(p_impl->levels[d + 1]))[v];
        s_cofaces.push_back(p_impl->handle_to_cell(*face_h));
    }
    return s_cofaces;
}

chain simplicial_complex::new_dense_chain(int d) {
    std::vector<double> v;
    if (d <= dimension())
        v = std::vector<double>(get_level_size(d), 0);
    else
        v = std::vector<double>(0, 0);
    return chain(d, v);
}

chain simplicial_complex::new_sparse_chain(int d) {
    vector_t v;
    if (d <= dimension())
        v = vector_t(get_level_size(d));
    else
        v = vector_t(0);
    return chain(d, v);
}

chain simplicial_complex::new_chain(int d) { return new_dense_chain(d); }

void simplicial_complex::calculate_hasse() {
    p_impl->has_hasse = true;
    p_impl->incidence = hasse_diag(*this);
}

// takes a dimension d
// returns a d-chain where the i-th entry is the volume of the i-th d-cell
chain simplicial_complex::volume_chain(int d) {
    chain volume_chain = new_dense_chain(d);
    for (cell_t cell : get_level(d)) {
        double volume;
        if (d == 0)
            volume = 1;
        else {
            point_t vec;
            point_t b_point = p_impl->points[cell[0]];
            vector<vector<double>> cell_edges;

            for (int i = 1; i < cell.size(); i++) {
                point_t i_point = p_impl->points[cell[i]];
                point_t n_edge;
                for (int j = 0; j < b_point.size(); j++) {
                    double v = i_point[j] - b_point[j];
                    n_edge.push_back(v);
                }
                cell_edges.push_back(n_edge);
            }

            ext_vector v(cell_edges[0]);
            for (int i = 1; i < cell_edges.size(); i++) {
                ext_vector w(cell_edges[i]);
                v = v ^ w;
            }
            volume = v.norm() / factorial(v.algebraic_dimension());
        }
        size_t cell_ind = cell_to_index(cell);
        volume_chain[cell_ind] = volume;
    }
    return volume_chain;
}
};
