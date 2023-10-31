/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : branch_bipartite_core_maintenance.h
* @details    : a fast bipartite core maintenance method
* @version    : 1.0
* @date       : 2022/01/30
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core_left_store_index.h"
#include "bipartite_core/bipartite_core_right_store_index.h"
#include "bipartite_core/bipartite_core_branch_store_index.h"
#include "bipartite_core/bipartite_core_order_index.h"
#include "bipartite_core/bipartite_core_rem_degree_index.h"
#include "bipartite_core/bipartite_core_degree_index.h"

namespace scnu{
    class branch_bipartite_core_maintenance {
    public:
        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map);

        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                         const shared_ptr<thread_pool> &pool);

        static void init(const shared_ptr<abstract_bipartite_graph> &B,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                         const shared_ptr<thread_pool>& pool);

        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map);

        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<bipartite_core_order_index> &core_order_index,
                           const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                           const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                           const shared_ptr<thread_pool> &pool);

        static void insert2(const shared_ptr<abstract_bipartite_graph> &B,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                            const shared_ptr<bipartite_core_order_index> &core_order_index,
                            const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                            const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                            const shared_ptr<thread_pool> &pool);


        static void insert(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<thread_pool> &pool);


        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map);

        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<bipartite_core_order_index> &core_order_index,
                           const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                           const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                           const shared_ptr<thread_pool> &pool);

        static void remove2(const shared_ptr<abstract_bipartite_graph> &B,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                            const shared_ptr<bipartite_core_order_index> &core_order_index,
                            const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                            const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                            const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<abstract_bipartite_graph> &B,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                           const shared_ptr<thread_pool> &pool);
    private:


        static uint32_t compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &B,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_vertex_index_map,
                                                        uint32_t l,
                                                        uint32_t i,
                                                        uint32_t j);

        static uint32_t compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &B,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_vertex_index_map,
                                                         uint32_t r,
                                                         uint32_t i,
                                                         uint32_t j);

        static uint32_t compute_left_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                        uint32_t l,
                                                        uint32_t i,
                                                        uint32_t j);

        static uint32_t compute_right_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                         uint32_t r,
                                                         uint32_t i,
                                                         uint32_t j);

        static uint32_t compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &B,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_vertex_index_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                        uint32_t l,
                                                        uint32_t i,
                                                        uint32_t j);

        static uint32_t compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &B,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_vertex_index_map,
                                                        const shared_ptr<unordered_set<uint32_t>> &removed_r_set,
                                                        uint32_t l,
                                                        uint32_t i,
                                                        uint32_t j);

        static uint32_t compute_left_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_vertex_index_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                        uint32_t l,
                                                        uint32_t i,
                                                        uint32_t j);


        static uint32_t compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &B,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_vertex_index_map,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                         uint32_t r,
                                                         uint32_t i,
                                                         uint32_t j);

        static uint32_t compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &B,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_vertex_index_map,
                                                         const shared_ptr<unordered_set<uint32_t>>& removed_l_set,
                                                         uint32_t r,
                                                         uint32_t i,
                                                         uint32_t j);

        static uint32_t compute_right_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_vertex_index_map,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                         uint32_t r,
                                                         uint32_t i,
                                                         uint32_t j);

        static void left_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                                         uint32_t i,
                                         uint32_t k,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map);

        static void left_candidate_graph(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                         const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                         const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_order_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                         uint32_t i,
                                         uint32_t k);

        static void left_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                uint32_t i,
                uint32_t k,
                const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map);


        static void middle_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                                           uint32_t k,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map);

        static void middle_candidate_graph(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                           const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                           const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_order_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                           uint32_t k);

        static void candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                    uint32_t k,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map);


        static void candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_order_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_degree_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                    uint32_t k,
                                    const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                    const shared_ptr<unordered_set<uint32_t>> &inserted_r_set);

        static void middle_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                uint32_t k,
                                           const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                                           const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map);


        static void right_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                                          uint32_t k,
                                          uint32_t j,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map);

        static void right_candidate_graph(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                          const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                          const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                          uint32_t k,
                                          uint32_t j);

        static void right_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                uint32_t k,
                uint32_t j,
                                          const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                                          const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map);


        static void left_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map,
                                                    uint32_t k);

        static void left_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<mutex>& global_left_mutex,
                                                    const shared_ptr<mutex>& global_right_mutex,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                    uint32_t k);

        static void left_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                    const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                    const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                    uint32_t k);

        static void left_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                    const shared_ptr<mutex>& global_left_mutex,
                                                    const shared_ptr<mutex>& global_right_mutex,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map,
                                                    uint32_t k);

        static void middle_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                      const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                                                      const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map,
                                                      const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map,
                                                       uint32_t k);

        static void middle_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                      const shared_ptr<mutex>& global_left_mutex,
                                                      const shared_ptr<mutex>& global_right_mutex,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                      uint32_t k,
                                                      const shared_ptr<thread_pool> &pool);

        static void middle_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                      const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                      const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                      const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                      const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                      const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                      uint32_t max_k,
                                                      const shared_ptr<thread_pool> &pool);

        static void middle_partial_core_decomposition2(const shared_ptr<abstract_bipartite_graph> &B,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                       const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                       const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                       const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                       const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                       const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                       uint32_t max_k,
                                                       const shared_ptr<thread_pool> &pool);

        static void middle_partial_core_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                const shared_ptr<mutex> &global_left_mutex,
                const shared_ptr<mutex> &global_right_mutex,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                uint32_t k,
                                                      const shared_ptr<thread_pool> &pool);

        static void right_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                      const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                                                      const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map,
                                                      uint32_t k);

        static void right_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                     const shared_ptr<mutex>& global_left_mutex,
                                                     const shared_ptr<mutex>& global_right_mutex,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                     uint32_t k);

        static void right_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                     const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                     const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                     uint32_t k);

        static void right_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                     const shared_ptr<mutex>& global_left_mutex,
                                                     const shared_ptr<mutex>& global_right_mutex,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                     uint32_t k);

        static void left_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                              const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                              const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                              uint32_t k,
                                              uint32_t j,
                                              const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                              const shared_ptr<unordered_set<uint32_t>> &removed_r_set);

        static void left_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                              const shared_ptr<vector<uint32_t>> &previous_removed_vector,
                                              uint32_t i,
                                              uint32_t k,
                                              const shared_ptr<vector<uint32_t>> &removed_vector);

        static void left_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                              const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                              const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                              uint32_t i,
                                              uint32_t k,
                                              const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                              const shared_ptr<unordered_set<uint32_t>> &removed_r_set);

        static void middle_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                uint32_t k,
                                                const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                const shared_ptr<unordered_set<uint32_t>> &removed_r_set);

        static void middle_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<vector<uint32_t>> &previous_removed_vector,
                                                uint32_t k,
                                                const shared_ptr<vector<uint32_t>> &removed_vector);

        static void middle_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                uint32_t k,
                                                const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                const shared_ptr<unordered_set<uint32_t>> &removed_r_set);

        static void right_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                               const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                               const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                               uint32_t i,
                                               uint32_t k,
                                               const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                               const shared_ptr<unordered_set<uint32_t>> &removed_r_set);

        static void right_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                               const shared_ptr<vector<uint32_t>> &previous_removed_vector,
                                               uint32_t k,
                                               uint32_t j,
                                               const shared_ptr<vector<uint32_t>> &removed_vector);

        static void right_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                               const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                               const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                               uint32_t i,
                                               uint32_t k,
                                               const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                               const shared_ptr<unordered_set<uint32_t>> &removed_r_set);

        static void remove_left_vertex(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_l_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_r_map,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                       uint32_t l,
                                       uint32_t i,
                                       uint32_t j);

        static void remove_left_vertex(const shared_ptr<abstract_bipartite_graph> &B,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_l_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_r_map,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_l_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_r_map,
                                       uint32_t l,
                                       uint32_t i,
                                       uint32_t j);

        static void remove_left_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_degree_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_degree_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                       uint32_t l,
                                       uint32_t i,
                                       uint32_t j);

        static void remove_left_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_degree_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_degree_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                       const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_l_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_r_map,
                                       uint32_t l,
                                       uint32_t i,
                                       uint32_t j);

        static void remove_right_vertex(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_r_map,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                        uint32_t r,
                                        uint32_t i,
                                        uint32_t j);

        static void remove_right_vertex(const shared_ptr<abstract_bipartite_graph> &B,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_r_map,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_r_map,
                                        uint32_t r,
                                        uint32_t i,
                                        uint32_t j);


        static void remove_right_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                        uint32_t r,
                                        uint32_t i,
                                        uint32_t j);

        static void remove_right_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                        const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_l_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_r_map,
                                        uint32_t r,
                                        uint32_t i,
                                        uint32_t j);

        static void remove_vertex(const shared_ptr<abstract_bipartite_graph> &B,
                                  uint32_t v,
                                  const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                  const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                  const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                  const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                  const shared_ptr<extend_list<int, uint32_t>> &order_list,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>>& k_rem_map,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &ext,
                                  const shared_ptr<map<int, uint32_t>>& affected_vertex_map,
                                  uint32_t i,
                                  uint32_t j);

        static void update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   uint32_t max_i,
                                                   uint32_t max_j);

        static void update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<mutex> &global_left_mutex,
                                                   const shared_ptr<mutex> &global_right_mutex,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   uint32_t max_i,
                                                   uint32_t max_j,
                                                   const shared_ptr<thread_pool> &pool);

        static void update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map);

        static void update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<thread_pool> &pool);


        static void update_index_for_insertion(const shared_ptr<abstract_bipartite_graph> &B,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                               const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                               const shared_ptr<unordered_set<uint32_t>> &inserted_r_set);

        static void update_index_for_removal(const shared_ptr<abstract_bipartite_graph> &B,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                             const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                             const shared_ptr<unordered_set<uint32_t>> &removed_r_set);

        static void update_left_core_degree_map(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &core_degree_map,
                                                uint32_t i,
                                                uint32_t k);

        static void update_middle_core_degree_map(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_left_store_index>>> &left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_right_store_index>>> &right_index_map,
                                                  const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                  const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &middle_core_degree_map,
                                                  uint32_t k);

        static void update_right_core_degree_map(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_left_store_index>>> &left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_right_store_index>>> &right_index_map,
                                                 const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                 const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &core_degree_map,
                                                 uint32_t k, uint32_t j);
    };
}

