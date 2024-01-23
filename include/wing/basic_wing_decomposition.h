/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : base_wing_decomposition.h
* @brief      : Basic decomposition methods
* @version    : 1.0
* @date       : 2020/9/7
******************************************************************************************************************/

#pragma once
#include "wing_utility.h"
#include "bipartite_edge_index.h"

namespace scnu
{
    /**
     * @brief basic decomposition method
     * @cite Zou Z. Bitruss decomposition of bipartite graphs[C] //International Conference on Database Systems for Advanced
     * Applications. Springer, Cham, 2016: 218-233.
     */
    class basic_wing_decomposition {
    public:

        static void init(const shared_ptr<abstract_bipartite_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         const shared_ptr<thread_pool> &pool);

        static void init(const shared_ptr<abstract_bipartite_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                         const shared_ptr<thread_pool> &pool);

        static void init(const shared_ptr<abstract_bipartite_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                         const shared_ptr<thread_pool> &pool);

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                  const shared_ptr<thread_pool>& pool);

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_support_map,
                                  const shared_ptr<thread_pool> &pool);

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                  const shared_ptr<thread_pool>& pool);


        static void edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map);

        static void edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map);

        static void edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                             const shared_ptr<thread_pool> &pool);


        static void vertex_priority_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_priority_map);


        static void set_edge_rank(const shared_ptr<unordered_set<shared_ptr<scnu::abstract_bipartite_edge>>> &edge_set,
                                  const shared_ptr<unordered_map<shared_ptr<scnu::abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                  const shared_ptr<uint32_t> &rank_id);
    };
}




