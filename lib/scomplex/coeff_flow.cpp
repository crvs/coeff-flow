#include <tuple>
#include <queue>
#include <scomplex/types.hpp>

namespace gsimp {

typedef std::tuple<cell_t, cell_t, double> q_elem_t;
typedef std::queue<q_element_t> queue_t;

class out_of_context : std::exception {};
class no_bounding_chain : std::exception {};

int chain_dim(chain_t p) { return std::get<0>(p); }
vector_t chain_rep(chain_t p) { return std::get<1>(p); }

chain_t coeff_flow(simplicial_complex& s_comp,  //
                   chain_t p,                   //
                   simplex_t sigma_0,           //
                   double c_0) {                //

    if (chain_dim(p) != s_comp.dim() - 1) throw out_of_context();

    std::vector<bool> seen_sigma(s_comp.level_size(s_comp.dim()), false);
    std::vector<bool> seen_tau(s_comp.level_size(s_comp.dim() - 1), false);
    vector_t b_chain;
    queue_t queue;

    cell_t sigma, tau;
    double c;
    size_t sigma_i, tau_i;
    cell_t while (not queue.empty()) {
        std::tie(sigma, tau, c) = queue.pop();
        sigma_i = s_comp.get_index_of_simplex(sigma);
        tau_i = s_comp.get_index_of_simplex(tau);
        if (seen_sigma.at(sigma_i)) {
            if (b_chain.coeffRef(sigma_i) != c) throw no_bounding_chain;
        } else {
            seen_tau.at(tau_i) = true;
            seen_sigma.at(sigma_i) = true;
            b_chain.coeffRef(sigma_i) = c;
            std::vector<cell_t> tau_cofaces;
            for (auto coface : s_comp.cofaces(tau)) {
                size_t coface_i = s_comp.get_index_of_simplex(coface);
                if (coface_i != sigma_i) tau_cofaces.push_back(coface);
            }
            // need cofaces in index form
        }
    }
}
};  // namespace gsimp
