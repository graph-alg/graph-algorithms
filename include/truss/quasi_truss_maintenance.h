/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : quasi_truss_maintenance.h
* @brief      : a parallel truss maintenance
* @version    : 1.0
* @date       : 2023/04/21
******************************************************************************************************************/
#pragma once

#include "truss/truss_utility.h"

namespace scnu {
    class quasi_truss_maintenance {
        static void init(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<thread_pool> &pool);

        static void init(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                         const shared_ptr<thread_pool> &pool);


        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                           const shared_ptr<uint32_t> &previous_k_max,
                           const shared_ptr<thread_pool> &pool);

        static void insert2(const shared_ptr<abstract_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                            const shared_ptr<uint32_t> &previous_k_max,
                            const shared_ptr<thread_pool> &pool);

        static void insert3(const shared_ptr<abstract_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_wing_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &WS,
                            const shared_ptr<uint32_t> &previous_k_max,
                            const shared_ptr<thread_pool> &pool);


        static void remove(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                           const shared_ptr<uint32_t> &previous_k_max,
                           const shared_ptr<thread_pool> &pool);

        static void remove2(const shared_ptr<abstract_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                            const shared_ptr<uint32_t> &previous_k_max,
                            const shared_ptr<thread_pool> &pool);

        static void remove3(const shared_ptr<abstract_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &WS,
                            const shared_ptr<uint32_t> &previous_k_max,
                            const shared_ptr<thread_pool> &pool);

    private:
        static void candidate_graph_finding(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                            uint32_t k,
                                            const shared_ptr<thread_pool> &pool);

        static void candidate_graph_finding2(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static uint32_t edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<abstract_edge> &e,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                 uint32_t k);


        static uint32_t edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<abstract_edge> &e,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_wing_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                 uint32_t k);

        static void edge_support_computation(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                             const shared_ptr<thread_pool> &pool);


        static uint32_t partial_truss_decomposition(const shared_ptr<abstract_graph> &G,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                    uint32_t k,
                                                    const shared_ptr<thread_pool> &pool);

        static uint32_t partial_truss_decomposition2(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     uint32_t k,
                                                     const shared_ptr<thread_pool> &pool);


        static void remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);


        static void set_edge_rank(const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<uint32_t> &rank_id);

        static void update_single_truss(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                        uint32_t k,
                                        const shared_ptr<thread_pool> &pool);


        static void update_edge_truss_support(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                              uint32_t k,
                                              const shared_ptr<thread_pool> &pool);

        static void update_single_truss2(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                         uint32_t k,
                                         const shared_ptr<thread_pool> &pool);

    };
}

