/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : branch_bipartite_core_decomposition.h
* @details    : an branch bipartite core decomposition
* @version    : 1.0
* @date       : 2022/01/02
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core_left_store_index.h"
#include "bipartite_core/bipartite_core_right_store_index.h"
#include "bipartite_core/bipartite_core_branch_store_index.h"
#include "bipartite_core/bipartite_core_degree_index.h"
#include "bipartite_core/bipartite_core_order_index.h"
#include "bipartite_core/bipartite_core_rem_degree_index.h"

namespace scnu{
    class branch_bipartite_core_decomposition {
    public:

        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                         const shared_ptr<::thread_pool> &pool);

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map);


        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                  const shared_ptr<thread_pool> &pool);

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &right_index_map,
                                  const shared_ptr<thread_pool> &pool);

        static uint32_t decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                  const shared_ptr<bipartite_core_order_index> &core_order_index,
                                  const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                  const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                  const shared_ptr<thread_pool> &pool);

        static uint32_t decompose2(const shared_ptr<abstract_bipartite_graph> &B,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                   const shared_ptr<bipartite_core_order_index> &core_order_index,
                                   const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                   const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                   const shared_ptr<thread_pool> &pool);


        static void find_core(const shared_ptr<abstract_bipartite_graph> &B,
                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_degree_map,
                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_degree_map,
                              uint32_t k);

        static void find_core(const shared_ptr<abstract_bipartite_graph> &B,
                              const shared_ptr<unordered_set<uint32_t>> &l_set,
                              const shared_ptr<unordered_set<uint32_t>> &r_set,
                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>& middle_order_map,
                              const shared_ptr<unordered_map<uint32_t, uint32_t>>& middle_rem_degree_map,
                              const shared_ptr<unordered_map<uint32_t, uint32_t>>& middle_core_degree_map,
                              uint32_t k);


        static uint32_t find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                       uint32_t k);



        static uint32_t find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                        uint32_t k);


        static uint32_t find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                       uint32_t k);

        static uint32_t find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                       const shared_ptr<unordered_set<uint32_t>> &l_set,
                                       const shared_ptr<unordered_set<uint32_t>> &r_set,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_oder_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_degree_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                       uint32_t k);

        static uint32_t find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &left_index_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &right_index_map,
                                       uint32_t k);

        static uint32_t find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                        uint32_t k);

        static uint32_t find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                        const shared_ptr<unordered_set<uint32_t>> &l_set,
                                        const shared_ptr<unordered_set<uint32_t>> &r_set,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_oder_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                        uint32_t k);

        static uint32_t find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &left_index_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_branch_store_index>>> &right_index_map,
                                        uint32_t k);
    };
}
