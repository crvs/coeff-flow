#include <tuple>
#include <queue>
#include <utility>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>

namespace gsimp {

typedef std::tuple<cell_t, cell_t, double> q_elem_t;
typedef std::queue<q_elem_t> queue_t;

class out_of_context : std::exception {};
class no_bounding_chain : std::exception {};

chain_t coeff_flow(simplicial_complex& s_comp,  //
                   chain_t p,                   //
                   cell_t sigma_0,              //
                   double c_0) {                //

    if (chain_dim(p) != s_comp.dimension() - 1) throw out_of_context();
    // 01
    chain_t c_chain(s_comp.dimension(),
                    vector_t(s_comp.get_level_size(s_comp.dimension())));
    std::vector<bool> seen_sigma(s_comp.get_level_size(s_comp.dimension()),
                                 false);
    std::vector<bool> seen_tau(s_comp.get_level_size(s_comp.dimension() - 1),
                               false);

    // 04
    size_t tau_i, sigma_i(s_comp.cell_to_index(sigma_0));
    chain_val(c_chain, sigma_i) = c_0;
    seen_sigma.at(sigma_i) = true;

    // 05
    queue_t queue;
    for (cell_t tau : s_comp.cell_boundary(sigma_0))
        queue.emplace(sigma_0, tau, c_0);

    // 08
    while (not queue.empty()) {
        cell_t sigma, tau;
        double c;
        std::tie(sigma, tau, c) = queue.front();
        queue.pop();
        sigma_i = s_comp.cell_to_index(sigma);
        tau_i = s_comp.cell_to_index(tau);

        // 10
        if (seen_sigma.at(sigma_i) && chain_val(c_chain, sigma_i) != c)
            throw no_bounding_chain();
        else if (not seen_sigma.at(sigma_i)) {
            // 14
            seen_tau.at(tau_i) = true;
            seen_sigma.at(sigma_i) = true;
            chain_val(c_chain, sigma_i) = c;

            // 15
            std::vector<std::pair<int, cell_t>> tau_cofaces;
            for (auto coface : s_comp.get_cof_and_ind(tau)) {
                size_t coface_i = s_comp.cell_to_index(std::get<1>(coface));
                if (coface_i != sigma_i) tau_cofaces.push_back(coface);
            }

            // 16
            double predicted_bdry =
                s_comp.boundary_inclusion_index(tau, sigma) * c;
            if (tau_cofaces.size() == 0 &&              //
                predicted_bdry != chain_val(p, tau_i))  //
                throw no_bounding_chain();              // 18
            else if (tau_cofaces.size() != 0) {
                int sigma_p_ind;
                cell_t sigma_p;
                std::tie(sigma_p_ind, sigma_p) = tau_cofaces[0];  // 20
                double c_p =
                    (sigma_p_ind *
                     (chain_val(p, tau_i) -
                      s_comp.boundary_inclusion_index(tau, sigma) * c));
                for (cell_t tau_p : s_comp.cell_boundary(sigma_p)) {
                    if (not seen_tau.at(s_comp.cell_to_index(tau_p))) {
                        queue.emplace(sigma_p, tau_p, c_p);
                        seen_tau.at(tau_i) = true;
                    }
                }
            }
        }
    }
    return c_chain;
}
};  // namespace gsimp
