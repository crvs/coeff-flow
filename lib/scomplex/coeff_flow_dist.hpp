#pragma once

#include <tuple>
#include <iostream>
#include <queue>
#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>

// multithreading
#include <ctpl.h>
#include <mutex>

namespace gsimp {

typedef std::tuple<cell_t, cell_t, double> q_elem_t;
typedef std::queue<q_elem_t> queue_t;

class out_of_context : std::exception {};
class no_bounding_chain : std::exception {};

chain_t coeff_flow(simplicial_complex& s_comp,  //
                   chain_t p,                   //
                   cell_t sigma_0,              //
                   double c_0,                  //
                   int num_cores) {             //

    if (chain_dim(p) != s_comp.dimension() - 1) throw out_of_context();
    // 01
    chain_t c_chain = s_comp.new_chain(s_comp.dimension());

    std::vector<bool> seen_sigma(
        (size_t)s_comp.get_level_size(s_comp.dimension()), false);

    std::vector<bool> seen_tau(
        (size_t)s_comp.get_level_size(s_comp.dimension() - 1), false);

    // 04
    size_t sigma_i(s_comp.cell_to_index(sigma_0));
    chain_val(c_chain, sigma_i) = c_0;
    seen_sigma[sigma_i] = true;

    // 05
    queue_t queue;
    for (cell_t tau : s_comp.cell_boundary(sigma_0)) {
        queue.emplace(sigma_0, tau, c_0);
        seen_tau[s_comp.cell_to_index(tau)] = true;
    }

    // threading facilities
    ctpl::thread_pool tp(num_cores);
    std::mutex mt0, mt1, mt2, mt3;
    // size_t num_taus = 0;

    // 08
    while (true) {
        std::cout << "woah";
        // num_taus++;
        // DEBUG OUTPUT
        /*
         * size_t ind = 0;
         * if (ind % 1000 == 0) {
         *     std::cout << ind << ": ";
         *     std::cout << " taus visited: " << seen_taus << "/" <<
         * seen_tau.size();
         *     std::cout << " sigmas visited: " << seen_sigmas << "/" <<
         * seen_sigma.size();
         *     std::cout << " queue size: " << queue.size() << "\n";
         *     std::cout.flush();
         * }
         * ind++;
         */

        if (tp.n_idle() > 0 && not(queue.empty()))
        {
        tp.push([&mt0,&mt1,&mt2,&queue,&s_comp,&p,&seen_sigma,&seen_tau,&c_chain](int id){
                    cell_t sigma;
                    cell_t tau;
                    double c;

                    mt2.lock();
                    std::tie(sigma, tau, c) = queue.front();
                    queue.pop();
                    mt2.unlock();
                    // get the indices of the respective cells
                    size_t sigma_i = s_comp.cell_to_index(sigma);
                    size_t tau_i = s_comp.cell_to_index(tau);

                    // check for local incoherence
                    mt0.lock();
                    if (seen_sigma[sigma_i] && chain_val(c_chain, sigma_i) != c)
                        throw no_bounding_chain();
                    else if (not(seen_sigma[sigma_i])) {
                        seen_sigma[sigma_i] = true;
                        chain_val(c_chain, sigma_i) = c;
                    }
                    mt0.unlock();

                    std::vector<std::pair<int, cell_t>> tau_cofaces;
                    for (auto coface : s_comp.get_cof_and_ind(tau)) {
                        size_t coface_i = s_comp.cell_to_index(std::get<1>(coface));
                        if (coface_i != sigma_i) tau_cofaces.push_back(coface);
                    }
                    double pi = p.second.coeffRef(tau_i);
                    double predicted_bdry = s_comp.boundary_inclusion_index(tau, sigma) * c;
                    if (tau_cofaces.size() == 0 && predicted_bdry != pi)
                        throw no_bounding_chain();
                    else if (tau_cofaces.size() != 0) {
                        int sigma_p_ind;
                        cell_t sigma_p;
                        std::tie(sigma_p_ind, sigma_p) = tau_cofaces[0];
                        double c_p = sigma_p_ind * (pi - s_comp.boundary_inclusion_index(tau, sigma) * c);
                        mt1.lock();
                        for (cell_t tau_p : s_comp.cell_boundary(sigma_p)) {
                            if (not seen_tau[s_comp.cell_to_index(tau_p)]) {
                                seen_tau[s_comp.cell_to_index(tau_p)] = true;
                                queue.emplace(sigma_p, tau_p, c_p);
                            }
                        }
                        mt1.unlock();
                    }
                });

        }
        else if (tp.size() == 0 && queue.empty()) break;
        }
    // std::cout << "number of taus: "<< num_taus << '\n';
    return c_chain;
}

/*
chain_t coeff_flow_embedded(simplicial_complex& s_comp, chain_t p) {
    if (chain_dim(p) != s_comp.dimension() - 1) throw out_of_context();

    cell_t sigma;
    double c;
    auto level = s_comp.get_level(s_comp.dimension() - 1);
    for (auto tau : s_comp.get_level(s_comp.dimension() - 1)) {
        if (s_comp.get_cofaces(tau).size() < 2) {
            int sigma_ind;
            size_t tau_i = s_comp.cell_to_index(tau);
            std::tie(sigma_ind, sigma) = s_comp.get_cof_and_ind(tau)[0];
            c = sigma_ind * chain_val(p, tau_i);
            chain_t sol = coeff_flow(s_comp, p, sigma, c);
            return sol;
        }
    }
    throw out_of_context();
}
*/
};  // namespace gsimp
