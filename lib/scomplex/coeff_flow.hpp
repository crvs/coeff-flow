#pragma once

#include <tuple>
#include <iostream>
#include <queue>
#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include "types.hpp"
#include "simplicial_complex.hpp"

namespace gsimp {
using namespace std;

typedef tuple<cell_t, cell_t, double> q_elem_t;
typedef queue<q_elem_t> queue_t;

class out_of_context : exception {};
class no_bounding_chain : exception {};

chain_v coeff_flow(simplicial_complex s_comp,  //
                   chain_v p,                   //
                   const cell_t& sigma_0,              //
                   const double& c_0) {                //

    if (chain_dim(p) != s_comp.dimension() - 1) throw out_of_context();
    // 01
    vector<double> c_vec(s_comp.get_level_size(s_comp.dimension()),0);

    vector<bool> seen_sigma((size_t)s_comp.get_level_size(s_comp.dimension()),
                            false);
    vector<bool> seen_tau((size_t)s_comp.get_level_size(s_comp.dimension() - 1),
                          false);

    size_t seen_taus = 0;
    seen_sigma.at(s_comp.cell_to_index(sigma_0)) = true;
    c_vec.at(s_comp.cell_to_index(sigma_0)) = c_0;
    size_t seen_sigmas = 1;

    queue_t queue;
    for (cell_t tau : s_comp.cell_boundary(sigma_0)) {
        queue.emplace(sigma_0, tau, c_0);
    }

    while (not queue.empty()) {
        /* // useful debug output
         *
         * size_t ind = 0;
         * if (ind % 1000 == 0) {
         *     cout << ind << ": ";
         *     cout << " taus visited: " << seen_taus << "/" <<
         * seen_tau.size();
         *     cout << " sigmas visited: " << seen_sigmas << "/" <<
         * seen_sigma.size();
         *     cout << " queue size: " << queue.size() << "\n";
         *     cout.flush();
         * }
         * ind++;
         */

        // dequeue the first element
        cell_t sigma;
        cell_t tau;
        double c;
        tie(sigma, tau, c) = queue.front();
        queue.pop();

        // get the indices of the respective cells
        size_t sigma_i = s_comp.cell_to_index(sigma);
        size_t tau_i = s_comp.cell_to_index(tau);

        if (seen_sigma.at(sigma_i)) {
            // found local incoherence
            if (c_vec.at(sigma_i) != c) throw no_bounding_chain();
        } else {
            seen_sigma.at(sigma_i) = true;
            c_vec.at(sigma_i) = c;
            seen_sigmas++;
        }

        if (seen_tau.at(tau_i)) continue;
        seen_tau.at(tau_i) = true;
        seen_taus++;

        cell_t sigma_p;
        size_t sigma_p_i;
        bool is_boundary = true;
        for (auto coface : s_comp.get_cof_and_ind(tau)) {
            size_t coface_i = s_comp.cell_to_index(get<1>(coface));
            if (coface_i != sigma_i) {
                is_boundary = false;
                sigma_p_i = coface_i;
                sigma_p = s_comp.index_to_cell(s_comp.dimension(), sigma_p_i);
            }
        }

        double predicted_bdry = s_comp.boundary_inclusion_index(tau, sigma) * c;
        if (is_boundary) {
            // sigma is the only coface of tau
            // check that we get the same value on tau
            if (predicted_bdry != chain_val(p, tau_i))
                throw no_bounding_chain();
        } else {
            // there is another coface we now focus on it
            double c_p = s_comp.boundary_inclusion_index(tau, sigma_p) *
                         (chain_val(p, tau_i) -
                          s_comp.boundary_inclusion_index(tau, sigma) * c);
            for (cell_t tau_p : s_comp.cell_boundary(sigma_p)) {
                size_t tau_p_i = s_comp.cell_to_index(tau_p);
                if (not seen_tau.at(tau_p_i)) {
                    // each tau only needs to be processed once
                    queue.emplace(sigma_p, tau_p, c_p);
                }
        }
    }
    }

    /*
     * // useful live debugging info
     * cout << "sigmas seen: " << seen_sigmas << "/" << s_comp.get_level_size(s_comp.dimension());
     * cout << " taus seen: " << seen_taus << "/" << s_comp.get_level_size(s_comp.dimension() - 1);
     * cout << '\n';
     */

    chain_v c_chain(s_comp.dimension(),c_vec);
    return c_chain;
}

chain_v coeff_flow_embedded(simplicial_complex s_comp, chain_v p) {
    if (chain_dim(p) != s_comp.dimension() - 1) throw out_of_context();

    cell_t sigma;
    double c;
    for (auto tau : s_comp.get_level(s_comp.dimension() - 1)) {
        if (s_comp.get_cofaces(tau).size() < 2) {
            int sigma_ind;
            size_t tau_i = s_comp.cell_to_index(tau);
            tie(sigma_ind, sigma) = s_comp.get_cof_and_ind(tau).at(0);
            c = sigma_ind * chain_val(p, tau_i);
            chain_v sol = coeff_flow(s_comp, p, sigma, c);
            return sol;
        }
    }
    throw out_of_context();
}
};  // namespace gsimp
