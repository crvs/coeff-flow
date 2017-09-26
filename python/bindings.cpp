#include <Python.h>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python.hpp>
#include <boost/function.hpp>

#include <scomplex_ho/types.hpp>
#include <scomplex_ho/chains.hpp>
#include <scomplex_ho/simplicial_complex.hpp>
#include <scomplex_ho/path_snapper.hpp>
#include <scomplex_ho/coeff_flow.hpp>
#include <scomplex_ho/quotient.hpp>

#include <iostream>
#include <vector>

#include <thread>
#include <chrono>

/*
 * phase 1:
 * construction:
 * take points, simplices, and paths (from python)
 */

using namespace std;

namespace py = boost::python;

template <typename T>
inline vector<T> to_std_vector(const py::object& iterable) {
    return vector<T>(py::stl_input_iterator<T>(iterable),
                     py::stl_input_iterator<T>());
}

template <typename T>
inline py::list to_py_list(const vector<T>& v) {
    py::object get_iter = py::iterator<vector<T>>();
    py::object iter = get_iter(v);
    py::list l(iter);
    return l;
}

/*
 * template <typename T>
 * inline void print_vec(const vector<T>& v) {
 *     for ( auto c : v ) cout << c << " ";
 *     cout << "\n";
 * }
 */

using namespace gsimp;

class path_clusterer {
    typedef py::object point_o;
    typedef py::object cell_o;
    typedef py::object path_o;

    std::shared_ptr<simplicial_complex> s_comp;
    bool has_quotient;
    std::shared_ptr<quotient> quotient_ptr;
    std::shared_ptr<path_snapper> p_snap;

    vector<vector<point_t>> paths;
    vector<chain> chains;
    vector<vector<chain>> gen_chains;

    vector<vector<size_t>> hom_clusters;

    // quotient??
    void construct_clusters_and_chains() {
        vector<bool> marked(chains.size(), false);
        size_t l_unmarked = 0;
        while (!marked.at(l_unmarked)) {
            vector<chain> gen_cluster;
            vector<size_t> cluster;
            bool seen_unmarked = false;

            // mark current chain
            size_t current = l_unmarked;
            cluster.push_back(current);
            marked.at(current) = true;

            for (int i = current + 1; i < chains.size(); ++i) {
                // calculate bounding chain
                chain p = chains.at(current);
                p -= chains.at(i);

                try {
                    chain c = coeff_flow_embedded(*s_comp, p);
                    gen_cluster.push_back(quotient_ptr->unquotient_chain(c));
                    cluster.push_back(i);
                } catch (const no_bounding_chain& e) {
                    if (!seen_unmarked) {
                        cout << "new cluster!\n";
                        l_unmarked = i;
                        seen_unmarked = true;
                    }
                }
            }
            hom_clusters.push_back(cluster);
            gen_chains.push_back(gen_cluster);
        }
    }

    void builder_source(const py::object& points,  //
                        const py::object& cells) {
        // decode points list
        vector<point_o> points_o = to_std_vector<py::object>(points);
        vector<point_t> points_v;
        for (auto pt : points_o) points_v.push_back(to_std_vector<double>(pt));

        // decode cells list
        vector<cell_o> cells_o = to_std_vector<py::object>(cells);
        vector<cell_t> cells_v;
        for (auto cell : cells_o)
            cells_v.push_back(to_std_vector<size_t>(cell));

        // construct simplicial complex
        s_comp = std::shared_ptr<simplicial_complex>(
            new simplicial_complex(points_v, cells_v));

        // construct path snapper
        p_snap = std::shared_ptr<path_snapper>(new path_snapper(s_comp));
        has_quotient = false;
    }

   public:
    path_clusterer(const py::object& points,  //
                   const py::object& cells) {
        builder_source(points, cells);
    }

    path_clusterer(const py::object& points,  //
                   const py::object& cells,   //
                   const py::object& a_paths) {
        builder_source(points, cells);

        for (size_t i = 0; i < py::len(a_paths); ++i) {
            this->add_path(a_paths[i]);
        }

        // cunstruct the clusters
        construct_clusters_and_chains();
    }

    void add_path(const py::object& path) {
        vector<point_t> path_v;
        for (size_t i = 0; i < py::len(path); ++i)
            path_v.push_back(to_std_vector<double>(path[i]));
        chain path_chain = p_snap->snap_path_to_dense_chain(path_v);
        chain q_path_chain;
        if (has_quotient)
            q_path_chain = quotient_ptr->quotient_chain(path_chain);
        bool added = false;
        for (size_t i = 0; i < hom_clusters.size(); i++) {
            if (hom_clusters[i].size() > 0) try {
                    chain b_chain;
                    if (has_quotient) {
                        b_chain = coeff_flow_embedded(
                            *quotient_ptr->quotient_complex(),
                            chains[hom_clusters[i][0]] - q_path_chain);
                        b_chain = quotient_ptr->unquotient_chain(b_chain);
                    } else
                        b_chain = coeff_flow_embedded(
                            *s_comp, chains[hom_clusters[i][0]] - q_path_chain);
                    hom_clusters[i].push_back(paths.size());
                    gen_chains[i].push_back(b_chain);
                    added = true;
                    break;
                } catch (no_bounding_chain& e) {
                }
        }
        if (!added) {
            hom_clusters.push_back({paths.size()});
            gen_chains.push_back({});
        }

        /*
        chain p_chain = quotient_ptr->quotient_chain(path_chain);
        p_chain.to_dense();
        p_chain = quotient_ptr->unquotient_chain(p_chain);

        chain pp_chain = quotient_ptr->quotient_chain(p_chain);
        pp_chain.to_dense();

        pp_chain = quotient_ptr->unquotient_chain(pp_chain);
        for (size_t i = 0 ; i < p_chain.get_size() ; i++ ) {
            if ( (p_chain[i] != 0) || (path_chain[i] != 0) )
            {
                std::cout << i << ": ";
            cell_t p_cell = s_comp->index_to_cell(p_chain.dimension(),i);
            cell_t pp_cell = quotient_ptr->quotient_face(p_cell);
            for (auto j : p_cell)
                std::cout << j << " ";
            std::cout << "- ";
            for (auto j : pp_cell)
                std::cout << j << " ";
            std::cout << ": " << path_chain[i] << " " << p_chain[i] << '\n';
            }
        }
        */

        chains.push_back(path_chain);
        paths.push_back(path_v);
    };

    py::list get_quotient_faces(int lv) {
        std::vector<cell_t> faces =
            quotient_ptr->quotient_complex()->get_level(lv);
        std::vector<py::list> listed_faces;
        for (auto face : faces)
            listed_faces.push_back(to_py_list<size_t>(face));
        return to_py_list<py::list>(listed_faces);
    }

    /*
    path_clusterer(const py::object& points,  //
                   const py::object& cells,   //
                   const py::object& paths,   //
                   bool(quot_f)(point_t)) {
        path_clusterer(points, cells, paths);
        make_quotient(quot_f);
    }
    */

    chain get_chain(size_t i) { return chains[i]; }

    vector<size_t> total_g_chains() {
        vector<size_t> totals;
        for (auto c : gen_chains) totals.push_back(c.size());
        return totals;
    }

    boost::function<bool(point_t)> characteristic_function(py::object boxes) {
        vector<vector<vector<pair<bool, double>>>> n_boxes;
        for (py::ssize_t b_i = 0; b_i < py::len(boxes); b_i++) {
            std::cout << "check 1 : " << b_i << "\n";
            py::object box = boxes[b_i];
            std::cout << "check 2 : " << b_i << "\n";
            vector<vector<pair<bool, double>>> n_box;
            for (py::ssize_t c_i = 0; c_i < py::len(box); c_i++) {
                std::cout << "check 3 : " << b_i << c_i << "\n";
                py::object coord = box[c_i];
                std::cout << "check 4 : " << b_i << c_i << "\n";
                bool has_min, has_max;
                double c_min, c_max;
                try {
                    c_min = py::extract<double>(coord[0]);
                    std::cout << "check 5a : " << b_i << c_i << "\n";
                    has_min = true;
                } catch (...) {
                    std::cout << "check 5b : " << b_i << c_i << "\n";
                    has_min = false;
                }
                try {
                    c_max = py::extract<double>(coord[1]);
                    std::cout << "check 6a : " << b_i << c_i << "\n";
                    has_max = true;
                } catch (...) {
                    std::cout << "check 6b : " << b_i << c_i << "\n";
                    has_max = false;
                }
                vector<pair<bool, double>> n_coord;
                n_coord.push_back(pair<bool, double>(has_min, c_min));
                n_coord.push_back(pair<bool, double>(has_max, c_max));
                n_box.push_back(n_coord);
            }
            n_boxes.push_back(n_box);
        }
        std::cout << "finished making boxes\n";

        boost::function<bool(point_t)> char_fun([n_boxes](point_t pt) {
            bool is_in_boxes = false;
            for (auto box : n_boxes) {
                bool is_in_box = true;
                for (size_t c = 0; c < pt.size(); ++c) {
                    if (box[c][0].first)
                        is_in_box = is_in_box && (box[c][0].second < pt[c]);
                    if (box[c][1].first)
                        is_in_box = is_in_box && (box[c][1].second > pt[c]);
                    if (!is_in_box) break;
                }
                is_in_boxes = is_in_boxes || is_in_box;
                if (is_in_boxes) break;
            }
            return is_in_boxes;
        });
        return char_fun;
    }

    void make_quotient(boost::function<bool(point_t)> quot_f) {
        std::cout << "making quotient\n";
        std::cout.flush();
        quotient_ptr = std::shared_ptr<quotient>(new quotient(s_comp, quot_f));
        has_quotient = true;

        std::cout << "TESTING QUOTIENT:\n";
        std::cout.flush();
        for (int j = 0 ; j <= s_comp->dimension() ; j++ ) {
        for (auto cell : s_comp->get_level(j)) {
            std::cout << s_comp->cell_to_index(cell) << ": ";
            std::cout.flush();
            for (auto i : cell) std::cout << i << " ";

            auto q_cell = quotient_ptr->quotient_face(cell);
            std::cout << "|q| " << quotient_ptr->quotient_index(q_cell) << ": ";
            std::cout.flush();
            for (auto i : q_cell) std::cout << i << " ";

            auto b_cell = quotient_ptr->unquotient_face(q_cell);
            if ((q_cell.size() != 0) && (b_cell.size() != 0)) {
                std::cout << "|b| " << quotient_ptr->base_index(b_cell) << ": ";
                std::cout.flush();
                for (auto i : b_cell) std::cout << i << " ";
            }

            std::cout << "\n";
            std::cout.flush();
        }
        }
    }

    void make_quotient_aux(py::object boxes) {
        std::cout << "entered\n";
        boost::function<bool(point_t)> quot_f = characteristic_function(boxes);
        std::cout << "passed here\n";
        std::cout.flush();
        make_quotient(quot_f);
    }

    chain get_g_chain(size_t c, size_t i) { return gen_chains[c][i]; }

    size_t unquotient_pt(size_t pt) {
        cell_t pt_face = {pt};
        if (pt == 0 ) {
            pt_face.push_back(s_comp->get_points().size());
        } else {
            pt_face = quotient_ptr->unquotient_face(pt_face);
        }
        return pt_face[0];
    }

    py::list get_quotient_chain(size_t num) {
        chain num_chain = chains[num];
        num_chain = quotient_ptr->quotient_chain(num_chain);
        std::vector<py::list> my_chain;
        std::vector<cell_t> q_level = //
            quotient_ptr->quotient_complex()->get_level(num_chain.dimension());
        for ( size_t i = 0 ; i < num_chain.get_dense().size() ; i++ ) {
            if (std::abs(num_chain[i]) != 0 ) {
                cell_t q_cell = q_level[i];
                cell_t b_cell;
                for (auto j : q_cell)
                    b_cell.push_back(unquotient_pt(j));
                my_chain.push_back(to_py_list<size_t>(b_cell));
            }
        }
        return to_py_list<py::list>(my_chain);
    }


    py::list homology_clusters() {
        vector<py::list> clusters_ls;
        for (auto cluster : hom_clusters)
            clusters_ls.push_back(to_py_list<size_t>(cluster));
        py::list ret = to_py_list<py::list>(clusters_ls);
        return ret;
    }

    py::list get_face(int d, size_t i) {
        cell_t face = s_comp->index_to_cell(d, i);
        return to_py_list<size_t>(face);
    }

    py::list generating_chains() {
        vector<py::list> g_chains;
        for (auto c : gen_chains) g_chains.push_back(to_py_list<chain>(c));
        py::list ret = to_py_list<py::list>(g_chains);
        return ret;
    }

    void send_signal(int secs) {
        std::this_thread::sleep_for(std::chrono::seconds(secs));
    }
};

BOOST_PYTHON_MODULE(gsimp) {
    py::class_<vector<py::list>>("vector_l")                 //
        .def(py::vector_indexing_suite<vector<py::list>>())  //
        ;

    py::class_<vector<py::object>>("vector_o")                 //
        .def(py::vector_indexing_suite<vector<py::object>>())  //
        ;

    py::class_<point_t>("vector_d")                        //
        .def(py::vector_indexing_suite<vector<double>>())  //
        ;

    py::class_<vector<point_t>>("vector_vd")                       //
        .def(py::vector_indexing_suite<vector<vector<double>>>())  //
        ;

    py::class_<cell_t>("vector_i")                         //
        .def(py::vector_indexing_suite<vector<size_t>>())  //
        ;

    py::class_<vector<cell_t>>("vector_vi")                        //
        .def(py::vector_indexing_suite<vector<vector<size_t>>>())  //
        ;

    py::class_<chain>("Chain", py::init<int, point_t>())  //
        .def(py::init<int, point_t>())                    //
        .def("is_dense", &chain::is_dense)                //
        .def("is_sparse", &chain::is_sparse)              //
        .def("is_sparse", &chain::is_dense)               //
        .def("abs", &chain::abs)                          //
        .def("dimension", &chain::dimension)              //
        .def("get_dense", &chain::get_dense)              //
        .def("get_size", &chain::get_size)                //
        .def("get_sparse", &chain::get_sparse)            //
        .def("is_dense", &chain::is_dense)                //
        .def("is_sparse", &chain::is_sparse)              //
        .def("to_dense", &chain::to_dense)                //
        .def("to_sparse", &chain::to_sparse)              //
        .def(py::self * float())                          //
        .def(py::self + py::self)                         //
        .def(py::self - py::self)                         //
        .def(py::self += py::self)                        //
        .def(py::self -= py::self)                        //
        .def(py::self ^ py::self)                         //
        ;

    py::class_<path_clusterer>(                                          //
        "clusterer",                                                     //
        py::init<py::object, py::object>())                              //
        .def(py::init<py::object, py::object>())                         //
        .def("homology_clusters", &path_clusterer::homology_clusters)    //
        .def("generating_chains", &path_clusterer::generating_chains)    //
        .def("get_chain", &path_clusterer::get_chain)                    //
        .def("get_g_chain", &path_clusterer::get_g_chain)                //
        .def("total_g_chains", &path_clusterer::total_g_chains)          //
        .def("add_path", &path_clusterer::add_path)                      //
        .def("quotient", &path_clusterer::make_quotient_aux)             //
        .def("get_face", &path_clusterer::get_face)                      //
        .def("wait_some", &path_clusterer::send_signal)                  //
        .def("get_q_faces", &path_clusterer::get_quotient_faces)         //
        .def("unquotient_pt", &path_clusterer::unquotient_pt)            //
        .def("get_quotient_chain", &path_clusterer::get_quotient_chain)  //
        ;
}
