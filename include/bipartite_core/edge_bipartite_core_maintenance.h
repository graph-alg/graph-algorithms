/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : edge_bipartite_core_maintenance.h
* @details    : a head file for index bipartite core maintenance
* @version    : 1.0
* @date       : 2021/02/05
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core_utility.h"
#include "bipartite_core/bipartite_core.h"
#include "bipartite_core/bipartite_core_compare.h"

namespace scnu{
    class edge_bipartite_core_maintenance {
    public:

        static void batch_insert(const shared_ptr<abstract_bipartite_graph> &B,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                 const shared_ptr<uint32_t>& delta);

        static void batch_insert(const shared_ptr<abstract_bipartite_graph> &B,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                 const shared_ptr<uint32_t> &delta,
                                 const shared_ptr<thread_pool> &pool);

        static void batch_remove(const shared_ptr<abstract_bipartite_graph> &B,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                 const shared_ptr<uint32_t>& delta);

        static void batch_remove(const shared_ptr<abstract_bipartite_graph> &B,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                 const shared_ptr<uint32_t> &delta,
                                 const shared_ptr<thread_pool> &pool);

        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map);

        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                         const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<abstract_bipartite_edge> &e,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<uint32_t>& delta);

        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                           const shared_ptr<abstract_bipartite_edge> &e,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<uint32_t> &delta,
                           const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<abstract_bipartite_edge> &e,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<uint32_t>& delta);

        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                           const shared_ptr<abstract_bipartite_edge> &e,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<uint32_t> &delta,
                           const shared_ptr<thread_pool> &pool);

    private:
        static void add_left_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                        const shared_ptr<unordered_set<uint32_t>> &Tl,
                                        const shared_ptr<unordered_set<uint32_t>> &Tr,
                                        const shared_ptr<unordered_set<uint32_t>> &Cl,
                                        const shared_ptr<unordered_set<uint32_t>> &Cr,
                                        const shared_ptr<unordered_set<uint32_t>> &Sl,
                                        const shared_ptr<unordered_set<uint32_t>> &Sr,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_l,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_r,
                                        uint32_t l,
                                        uint32_t alpha,
                                        uint32_t beta);

        static void add_right_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                         const shared_ptr<unordered_set<uint32_t>> &Tl,
                                         const shared_ptr<unordered_set<uint32_t>> &Tr,
                                         const shared_ptr<unordered_set<uint32_t>> &Cl,
                                         const shared_ptr<unordered_set<uint32_t>> &Cr,
                                         const shared_ptr<unordered_set<uint32_t>> &Sl,
                                         const shared_ptr<unordered_set<uint32_t>> &Sr,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_l,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_r,
                                         uint32_t r,
                                         uint32_t alpha,
                                         uint32_t beta);

        static void batch_update_alpha_phi_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                      uint32_t alpha,
                                                      uint32_t phi_alpha);

        static void batch_update_alpha_phi_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                      uint32_t alpha,
                                                      uint32_t phi_alpha);

        static void batch_update_alpha_tau_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                      uint32_t alpha,
                                                      uint32_t tau_alpha);

        static void batch_update_alpha_tau_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                      uint32_t alpha,
                                                      uint32_t tau_alpha);

        static void batch_update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                    uint32_t phi_beta,
                                                    uint32_t beta);

        static void batch_update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                    uint32_t phi_beta,
                                                    uint32_t beta);

        static void batch_update_tau_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                    uint32_t tau_beta,
                                                    uint32_t beta);

        static void batch_update_tau_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                    uint32_t tau_beta,
                                                    uint32_t beta);

        static uint32_t find_core(const shared_ptr<abstract_bipartite_graph> &B);


        static uint32_t get_maximal_b_alpha(const shared_ptr<abstract_bipartite_graph> &B,
                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                     uint32_t l,
                                     uint32_t alpha);

        static uint32_t get_maximal_b_beta(const shared_ptr<abstract_bipartite_graph> &B,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                           uint32_t r,
                                           uint32_t beta);

        static void remove_left_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                           const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cl,
                                           const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cr,
                                           uint32_t l,
                                           uint32_t alpha,
                                           uint32_t beta);

        static void remove_right_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                            const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cl,
                                            const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cr,
                                            uint32_t r,
                                            uint32_t alpha,
                                            uint32_t beta);

        static void update_alpha_phi_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                const shared_ptr<abstract_bipartite_edge>& e,
                                                uint32_t alpha);

        static void update_alpha_phi_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                const shared_ptr<abstract_bipartite_edge> &e,
                                                uint32_t alpha);

        static void update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                              const shared_ptr<abstract_bipartite_edge>& e,
                                              uint32_t beta);


        static void update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                              const shared_ptr<abstract_bipartite_edge> &e,
                                              uint32_t beta);

        static void update_alpha_tau_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                const shared_ptr<abstract_bipartite_edge> &e,
                                                uint32_t alpha);

        static void update_alpha_tau_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                const shared_ptr<abstract_bipartite_edge> &e,
                                                uint32_t alpha);

        static void update_beta_tau_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                              const shared_ptr<abstract_bipartite_edge> &e,
                                              uint32_t beta);

        static void update_beta_tau_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                              const shared_ptr<abstract_bipartite_edge> &e,
                                              uint32_t beta);

    };
}



