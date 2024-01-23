/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : basic_truss_decomposition.h
* @brief      : common header files for k-truss
* @version    : 1.0
* @date       : 2022/6/17
******************************************************************************************************************/

#pragma once
#include "truss/truss_utility.h"

namespace scnu{
    class basic_truss_decomposition {
    public:
        static void init(const shared_ptr<scnu::abstract_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map);

        static void init(const shared_ptr<scnu::abstract_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map);

        static void init(const shared_ptr<scnu::abstract_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<thread_pool> &pool);

        static void init(const shared_ptr<scnu::abstract_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& rem,
                         const shared_ptr<thread_pool> &pool);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                  const shared_ptr<thread_pool>& pool);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                  const shared_ptr<thread_pool>& pool);

       static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_support_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>>& edge_mutex_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_rank_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_support_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& rem,
                                  const shared_ptr<thread_pool> &pool);

    private:

        static void assign_rank(const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                const shared_ptr<uint32_t>& base_rank);

        static void edge_support_computation(const shared_ptr<abstract_graph>& G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map);

        static void edge_support_computation(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<thread_pool> &pool);
    };
}


