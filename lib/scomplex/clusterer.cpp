#include <limits>

#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/chain_calc.hpp>
#include <scomplex/coeff_flow.hpp>
#include <scomplex/path_snapper.hpp>
#include <scomplex/clusterer.hpp>

namespace gsimp {
struct clusterer::impl {
    std::shared_ptr<simplicial_complex> s_comp;
    std::shared_ptr<path_snapper> p_snap;
    std::shared_ptr<bounding_chain> b_calc;

    std::vector<std::vector<point_t>> paths;
    std::vector<chain_t> path_chains;

    std::vector<int> homology_clusters;

    impl(simplicial_complex& sc) {
        s_comp =
            std::shared_ptr<simplicial_complex>(new simplicial_complex(sc));
        p_snap = std::shared_ptr<path_snapper>(new path_snapper(s_comp));
    }

    impl(path_snapper& psnap) {
        p_snap = std::shared_ptr<path_snapper>(new path_snapper(psnap));
        s_comp = p_snap->get_underlying_complex();
    }

    impl(std::shared_ptr<path_snapper> psnap) {
        p_snap = psnap;
        s_comp = p_snap->get_underlying_complex();
    }

    impl(std::shared_ptr<simplicial_complex> sc) {
        s_comp = sc;
        p_snap = std::shared_ptr<path_snapper>(new path_snapper(s_comp));
    }

    ~impl() {}

    bool have_homologies() { return homology_clusters.size() == paths.size(); }

    void calculate_homologies() {}

    void calculate_chains() {
        // get the bounding chain if it isn't there already
        b_calc = std::shared_ptr<bounding_chain>(new bounding_chain(s_comp));
        for (size_t i = 0; i < paths.size(); ++i) {
            for (size_t j = i + 1; j < paths.size(); ++j) {
                chain_t chain = subtract(path_chains[i], path_chains[i]);
                b_calc->get_bounding_chain(chain);
            }
        }
    }

    void add_path(std::vector<point_t> path) {
        paths.push_back(path);
        path_chains.push_back(p_snap->snap_path_to_chain(path));
    }
};

clusterer::clusterer(simplicial_complex& sc) {
    p_impl = std::shared_ptr<impl>(new impl(sc));
}
clusterer::clusterer(path_snapper& ps) {
    p_impl = std::shared_ptr<impl>(new impl(ps));
}

clusterer::clusterer(std::shared_ptr<simplicial_complex> sc) {
    p_impl = std::shared_ptr<impl>(new impl(sc));
}
clusterer::clusterer(std::shared_ptr<path_snapper> ps) {
    p_impl = std::shared_ptr<impl>(new impl(ps));
}

clusterer::clusterer(clusterer& other) { p_impl = other.p_impl; }
clusterer::~clusterer() {}

clusterer& clusterer::operator=(const clusterer& other) {
    p_impl = other.p_impl;
    return *this;
}

// adds a path to the clusterer object
void clusterer::add_path(std::vector<point_t> path) { p_impl->add_path(path); }

// returns a vector where the i-th entry, is the number of the homology
// cluster the i-th path is in
std::vector<int> clusterer::get_homology_clusters() {
    if (!p_impl->have_homologies()) p_impl->calculate_homologies();
    return p_impl->homology_clusters;
}

// returns the vector of pairwise distances between the paths in the
// clusterer object
std::vector<double> clusterer::get_dist_array() {}

std::vector<point_t> clusterer::get_path(size_t i) { return p_impl->paths[i]; }

chain_t clusterer::get_chain(size_t) { }

chain_t clusterer::get_bounding_chain(size_t, size_t) {}

};
