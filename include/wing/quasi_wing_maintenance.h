/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : quasi_wing_maintenance.h
* @version    : 1.1
* @date       : 2021/04/26
******************************************************************************************************************/

#pragma once

#include "wing_utility.h"
#include "BE_index.h"

namespace scnu {
    class quasi_wing_maintenance {
    public:

        static void previous_k_max(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                   uint32_t &k_max,
                                   const shared_ptr<thread_pool>& pool);




        static void init(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                         const shared_ptr<thread_pool>& pool);


        static void init(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_support_map,
                         const shared_ptr<thread_pool>& pool);


        static void insert(const shared_ptr<abstract_bipartite_graph> &G,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_support_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                           const shared_ptr<uint32_t>& previous_k_max,
                           const shared_ptr<thread_pool>& pool);

        static void insert2(const shared_ptr<abstract_bipartite_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_support_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                            const shared_ptr<uint32_t>& previous_k_max,
                            const shared_ptr<thread_pool>& pool);

        static void insert3(const shared_ptr<abstract_bipartite_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                            const shared_ptr<uint32_t> &previous_k_max,
                            const shared_ptr<thread_pool> &pool);

        static void insert4(const shared_ptr<abstract_bipartite_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_support_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                            const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                            const shared_ptr<uint32_t>& previous_k_max,
                            const shared_ptr<thread_pool>& pool);


        static void remove(const shared_ptr<abstract_bipartite_graph> &G,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                           const shared_ptr<uint32_t>& previous_k_max,
                           const shared_ptr<thread_pool>& pool);

        static void remove2(const shared_ptr<abstract_bipartite_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                            const shared_ptr<uint32_t>& previous_k_max,
                            const shared_ptr<thread_pool>& pool);


        static void remove3(const shared_ptr<abstract_bipartite_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                            const shared_ptr<uint32_t>& previous_k_max,
                            const shared_ptr<thread_pool>& pool);

        static void remove4(const shared_ptr<abstract_bipartite_graph> &G,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                            const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                            const shared_ptr<uint32_t>& previous_k_max,
                            const shared_ptr<thread_pool>& pool);

    private:


        static void candidate_graph_finding(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                            uint32_t k,
                                            const shared_ptr<thread_pool>& pool);

        static void candidate_graph_finding2(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                             uint32_t k,
                                             const shared_ptr<thread_pool>& pool);

        static void candidate_graph_finding3(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<BE_index>& bloom_index,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                             uint32_t k,
                                             const shared_ptr<thread_pool>& pool);

        static void candidate_graph_finding4(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                             const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                             uint32_t k,
                                             const shared_ptr<thread_pool>& pool);

        static uint32_t edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                 const shared_ptr<abstract_bipartite_edge> &e,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                 uint32_t k);


        static uint32_t edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                 const shared_ptr<abstract_bipartite_edge> &e,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                 uint32_t k);

        static void edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                             const shared_ptr<thread_pool> &pool);

        static shared_ptr<BE_index> index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                       const shared_ptr<thread_pool> &pool);

        static void left_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                            const shared_ptr<BE_index> &bloom_index,
                                            const shared_ptr<thread_pool>& pool);


        static uint32_t partial_wing_decomposition(const shared_ptr<abstract_bipartite_graph> &G,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                   uint32_t k,
                                                   const shared_ptr<thread_pool>& pool);

        static uint32_t partial_wing_decomposition2(const shared_ptr<abstract_bipartite_graph> &G,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                    uint32_t k,
                                                    const shared_ptr<thread_pool>& pool);

        static uint32_t partial_wing_decomposition3(const shared_ptr<abstract_bipartite_graph> &G,
                                                    const shared_ptr<BE_index> &bloom_index,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                    uint32_t k,
                                                    const shared_ptr<thread_pool> &pool);

        static uint32_t partial_wing_decomposition4(const shared_ptr<abstract_bipartite_graph> &G,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rank_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                    uint32_t k,
                                                    const shared_ptr<thread_pool>& pool);


        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                             const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);



        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &rectangle_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& next_edge_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &rectangle_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& next_edge_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                             const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<BE_index> &bloom_index,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<BE_index> &bloom_index,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& candidate_edge_support_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &rectangle_edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& next_edge_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void right_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_priority_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                             const shared_ptr<BE_index> &bloom_index,
                                             const shared_ptr<thread_pool>& pool);


        static void set_edge_rank(const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<uint32_t> &rank_id);

        static void update_single_wing(const shared_ptr<abstract_bipartite_graph> &G,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                       uint32_t k,
                                       const shared_ptr<thread_pool> &pool);

        static void update_single_wing2(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                        uint32_t k,
                                        const shared_ptr<thread_pool> &pool);

        static void update_single_wing3(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<BE_index> &bloom_index,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                        uint32_t k,
                                        const shared_ptr<thread_pool> &pool);

        static void update_single_wing4(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                        uint32_t k,
                                        const shared_ptr<thread_pool> &pool);



        static void update_edge_wing_support(const shared_ptr<abstract_bipartite_graph>& G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void update_edge_wing_support(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<BE_index> &bloom_index,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void update_edge_wing_support(const shared_ptr<abstract_bipartite_graph>& G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                             const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &candidate_edge_set,
                                             uint32_t k,
                                             const shared_ptr<thread_pool> &pool);

        static void wl_move(const shared_ptr<abstract_bipartite_edge> &e,
                            const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                            uint32_t source,
                            uint32_t destination);



        static void vertex_priority_computation(const shared_ptr<abstract_bipartite_graph>& G,
                                                const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_priority_map,
                                                const shared_ptr<thread_pool> &pool);

    };
}
