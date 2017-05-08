#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/chain_calc.hpp>
#include <scomplex/coeff_flow.hpp>
#include <scomplex/path_snapper.hpp>

namespace gsimp {
class clusterer {
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:
    clusterer(simplicial_complex&);
    clusterer(path_snapper&);
    clusterer(std::shared_ptr<simplicial_complex>);
    clusterer(std::shared_ptr<path_snapper>);

    clusterer(clusterer&);
    clusterer& operator=(const clusterer&);

    ~clusterer();

    // adds a path to the clusterer object
    void add_path(std::vector<point_t>);

    // returns a vector where the i-th entry, is the number of the homology
    // cluster the i-th path is in
    std::vector<int> get_homology_clusters();

    // returns the vector of pairwise distances between the paths in the
    // clusterer object
    std::vector<double> get_dist_array();

    std::vector<point_t> get_path(size_t);
    chain_t get_chain(size_t);
    chain_t get_bounding_chain(size_t, size_t);
};
};
