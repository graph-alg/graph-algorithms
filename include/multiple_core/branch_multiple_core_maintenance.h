/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : hierarchy_multiple_core_maintenance.h
* @brief      : a method for temporal core maintenance
* @version    : 1.0
* @date       : 2022/11/04
******************************************************************************************************************/

#pragma once

#include "multiple_core_utility.h"
#include "multiple_core_pair_map_index.h"

namespace scnu {
    class branch_multiple_core_maintenance {
    public:

        static void init(const shared_ptr<temporal_graph> &G,
                         const shared_ptr<unordered_map<uint32_t,
                                 shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_size_map,
                         const shared_ptr<thread_pool> &pool);


        static void init(const shared_ptr<temporal_graph> &G,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                         const shared_ptr<thread_pool> &pool);

        static void init(const shared_ptr<temporal_graph> &G,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                         const shared_ptr<unordered_map<uint32_t,
                                 shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_size_map,
                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                         const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t,
                                   shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                           const shared_ptr<thread_pool> &pool);

        static void insert(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t,
                                   shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                           const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_size_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                           const shared_ptr<thread_pool> &pool);

        static void remove(const shared_ptr<temporal_graph> &G,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_size_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                           const shared_ptr<thread_pool> &pool);

    private:

        static void assign(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                           uint32_t value,
                           const shared_ptr<thread_pool> &pool);

        static void insertion_assign(const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                     uint32_t k,
                                     uint32_t h,
                                     const shared_ptr<thread_pool>& pool);

        static void removal_assign(const shared_ptr<unordered_set<uint32_t>> &removed_set,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                   uint32_t k,
                                   uint32_t h,
                                   const shared_ptr<thread_pool>& pool);

        static void insertion_merge(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<thread_pool> &pool);


        static void removal_merge(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<thread_pool> &pool);

        static uint32_t compute_left_core_degree(const shared_ptr<temporal_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                 uint32_t u,
                                                 uint32_t k,
                                                 uint32_t h);

        static uint32_t compute_left_core_degree(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                uint32_t u,
                uint32_t k,
                uint32_t h);

        static uint32_t compute_middle_core_degree(const shared_ptr<temporal_graph> &G,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                   uint32_t u,
                                                   uint32_t k,
                                                   uint32_t h);

        static uint32_t compute_middle_core_degree(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                uint32_t u,
                uint32_t k,
                uint32_t h);

        static uint32_t compute_right_core_degree(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                  uint32_t u,
                                                  uint32_t k,
                                                  uint32_t h);

        static uint32_t compute_right_core_degree(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                uint32_t u,
                uint32_t k,
                uint32_t h);

        static uint32_t compute_core_degree(const shared_ptr<temporal_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                            uint32_t u,
                                            uint32_t k,
                                            uint32_t h);

        static uint32_t compute_core_degree(const shared_ptr<unordered_map<uint32_t,
                shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                            uint32_t u,
                                            uint32_t k,
                                            uint32_t h);

        static uint32_t find_max_delta(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map);

        static uint32_t find_max_delta(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                       const shared_ptr<thread_pool>& pool);

        static uint32_t find_max_delta(const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map);

        static uint32_t find_max_delta(const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                       const shared_ptr<thread_pool>& pool);
        static uint32_t find_max_k(const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                   uint32_t h);

        static uint32_t find_max_k(const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                   uint32_t h,
                                   const shared_ptr<thread_pool>& pool);


        static uint32_t find_max_h(const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                   uint32_t k);

        static uint32_t find_max_h(const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                   uint32_t k,
                                   const shared_ptr<thread_pool>& pool);

        static bool left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                         uint32_t k,
                                         uint32_t h);

        static bool left_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                uint32_t k,
                uint32_t h);

        static bool left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                         const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                         uint32_t k,
                                         uint32_t h,
                                         const shared_ptr<thread_pool> &pool);

        static bool left_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);

        static bool middle_candidate_graph(const shared_ptr<temporal_graph> &G,
                                           const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                           uint32_t k,
                                           uint32_t h);

        static bool middle_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                uint32_t k,
                uint32_t h);

        static bool middle_candidate_graph(const shared_ptr<temporal_graph> &G,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                           const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                           uint32_t k,
                                           uint32_t h,
                                           const shared_ptr<thread_pool> &pool);

        static bool middle_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);

        static bool right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                          const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                          uint32_t k,
                                          uint32_t h);

        static bool right_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                uint32_t k,
                uint32_t h);

        static bool right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                          const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                          uint32_t k,
                                          uint32_t h,
                                          const shared_ptr<thread_pool> &pool);

        static bool right_candidate_graph(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);

        static void left_decomposition(const shared_ptr<temporal_graph> &G,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                       uint32_t k,
                                       uint32_t h);

        static void left_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                uint32_t k,
                uint32_t h);


        static void left_decomposition(const shared_ptr<temporal_graph> &G,
                                       const shared_ptr<mutex> &global_mutex,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                       uint32_t k,
                                       uint32_t h);

        static void left_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<mutex> &global_mutex,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                uint32_t k,
                uint32_t h);

        static void left_decomposition(const shared_ptr<temporal_graph> &G,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                       uint32_t k,
                                       uint32_t h,
                                       const shared_ptr<thread_pool> &pool);

        static void left_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);


        static void middle_decomposition(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<mutex> &global_mutex,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                         uint32_t delta,
                                         const shared_ptr<thread_pool> &pool);

        static void middle_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<mutex> &global_mutex,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                uint32_t delta,
                const shared_ptr<thread_pool> &pool);

        static void middle_decomposition(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                         uint32_t delta,
                                         const shared_ptr<thread_pool> &pool);

        static void middle_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                uint32_t delta,
                const shared_ptr<thread_pool> &pool);


        static void right_decomposition(const shared_ptr<temporal_graph> &G,
                                        const shared_ptr<mutex> &global_mutex,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                        uint32_t k,
                                        uint32_t h);

        static void right_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<mutex> &global_mutex,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                uint32_t k,
                uint32_t h);

        static void right_decomposition(const shared_ptr<temporal_graph> &G,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                        uint32_t k,
                                        uint32_t h,
                                        const shared_ptr<thread_pool> &pool);

        static void right_decomposition(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);

        static void remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                uint32_t w,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                uint32_t k,
                                                uint32_t h);

        static void remove_unsatisfied_vertices(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                uint32_t w,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                uint32_t k,
                                                uint32_t h);

        static void remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                uint32_t k,
                                                uint32_t h);

        static void remove_unsatisfied_vertices(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                uint32_t k,
                uint32_t h);

        static void remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                uint32_t k,
                                                uint32_t h,
                                                const shared_ptr<thread_pool> &pool);

        static void remove_unsatisfied_vertices(const shared_ptr<unordered_map<uint32_t,
                shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                uint32_t k,
                                                uint32_t h,
                                                const shared_ptr<thread_pool> &pool);

        static shared_ptr<unordered_set<uint32_t>> left_removal_core(const shared_ptr<temporal_graph> &G,
                                                                     const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                     const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                                     uint32_t k,
                                                                     uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> left_removal_core(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                uint32_t k,
                uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> left_removal_core(const shared_ptr<temporal_graph> &G,
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                     const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                     const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                                     uint32_t k,
                                                                     uint32_t h,
                                                                     const shared_ptr<thread_pool> &pool);

        static shared_ptr<unordered_set<uint32_t>> left_removal_core(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);

        static shared_ptr<unordered_set<uint32_t>> middle_removal_core(const shared_ptr<temporal_graph> &G,
                                                                       const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                                       uint32_t k,
                                                                       uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> middle_removal_core(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                uint32_t k,
                uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> middle_removal_core(const shared_ptr<temporal_graph> &G,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                                       uint32_t k,
                                                                       uint32_t h,
                                                                       const shared_ptr<thread_pool> &pool);

        static shared_ptr<unordered_set<uint32_t>> middle_removal_core(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);

        static shared_ptr<unordered_set<uint32_t>> right_removal_core(const shared_ptr<temporal_graph> &G,
                                                                      const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                      const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                                      uint32_t k,
                                                                      uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> right_removal_core(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                uint32_t k,
                uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> right_removal_core(const shared_ptr<temporal_graph> &G,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                      const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                      const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                                      uint32_t k,
                                                                      uint32_t h,
                                                                      const shared_ptr<thread_pool> &pool);

        static shared_ptr<unordered_set<uint32_t>> right_removal_core(
                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                uint32_t k,
                uint32_t h,
                const shared_ptr<thread_pool> &pool);

        static void update_vertex_edge_index_for_insertion(const shared_ptr<scnu::temporal_graph> &G,
                                                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                           const shared_ptr<unordered_set<uint32_t>>& affected_set);

        static void update_vertex_edge_index_for_insertion(const shared_ptr<scnu::temporal_graph> &G,
                                                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                           const shared_ptr<unordered_set<uint32_t>>& affected_set,
                                                           const shared_ptr<thread_pool> &pool);

        static void update_vertex_edge_index_for_removal(const shared_ptr<scnu::temporal_graph> &G,
                                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                         const shared_ptr<unordered_set<uint32_t>>& affected_set,
                                                         const shared_ptr<unordered_set<uint32_t>>& isolated_vertex_set);

        static void update_vertex_edge_index_for_removal(const shared_ptr<scnu::temporal_graph> &G,
                                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                         const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                         const shared_ptr<unordered_set<uint32_t>> &isolated_vertex_set,
                                                         const shared_ptr<thread_pool> &pool);
    };
}