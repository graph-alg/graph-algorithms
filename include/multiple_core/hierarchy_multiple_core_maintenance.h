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

namespace scnu{
    class hierarchy_multiple_core_maintenance {
    public:

        static void init(const shared_ptr<temporal_graph>& G,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                         const shared_ptr<thread_pool>& pool);

        static void insert(const shared_ptr<temporal_graph>& G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map);

        static void insert(const shared_ptr<temporal_graph>& G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                           const shared_ptr<thread_pool>& pool);

        static void insert(const shared_ptr<temporal_graph>& G,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                           const shared_ptr<thread_pool>& pool);


        static void remove(const shared_ptr<temporal_graph>& G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map);

        static void remove(const shared_ptr<temporal_graph>& G,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                           const shared_ptr<thread_pool>& pool);

        static void remove(const shared_ptr<temporal_graph>& G,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>>& edge_set,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                           const shared_ptr<thread_pool>& pool);
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
                                                 uint32_t u,
                                                 uint32_t k,
                                                 uint32_t h);

        static uint32_t compute_right_core_degree(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  uint32_t u,
                                                  uint32_t k,
                                                  uint32_t h);

        static uint32_t compute_core_degree(const shared_ptr<temporal_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                            uint32_t u,
                                            uint32_t k,
                                            uint32_t h);

        static bool left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_map,
                                         uint32_t k,
                                         uint32_t h);

        static bool left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_map,
                                         uint32_t k,
                                         uint32_t h,
                                         const shared_ptr<thread_pool>& pool);

        static bool right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                          const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_map,
                                          uint32_t k,
                                          uint32_t h);

        static bool right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                          const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_map,
                                          uint32_t k,
                                          uint32_t h,
                                          const shared_ptr<thread_pool>& pool);

        static void find_left_quasi_core(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                         uint32_t k,
                                         uint32_t h);

        static void find_right_quasi_core(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                         uint32_t k,
                                         uint32_t h);

        static void find_left_quasi_core(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                         uint32_t k,
                                         uint32_t h,
                                         const shared_ptr<thread_pool>& pool);

        static void find_right_quasi_core(const shared_ptr<temporal_graph> &G,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                         const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                         uint32_t k,
                                         uint32_t h,
                                         const shared_ptr<thread_pool>& pool);

        static void find_quasi_cores(const shared_ptr<temporal_graph> &G,
                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                     const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_vertex_degree_map,
                                     uint32_t k,
                                     uint32_t h,
                                     const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>> &quasi_edge_map);

        static void find_quasi_cores(const shared_ptr<temporal_graph> &G,
                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                     const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_vertex_degree_map,
                                     uint32_t k,
                                     uint32_t h,
                                     const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>> &quasi_edge_map,
                                     const shared_ptr<thread_pool>& pool);


        static void remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                uint32_t w,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                uint32_t k,
                                                uint32_t h);

        static void remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                const shared_ptr<unordered_set<uint32_t>>& invalid_set,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                uint32_t k,
                                                uint32_t h,
                                                const shared_ptr<thread_pool>& pool);


        static bool update_single_core(const shared_ptr<temporal_graph> &G,
                                       const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                       const shared_ptr<unordered_set<uint32_t>>& removed_set,
                                       uint32_t k,
                                       uint32_t h);

        static bool  update_single_core(const shared_ptr<temporal_graph> &G,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                        const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                        const shared_ptr<unordered_set<uint32_t>>& removed_set,
                                        uint32_t k,
                                        uint32_t h,
                                        const shared_ptr<thread_pool>& pool);

    };
}



