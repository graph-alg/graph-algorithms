/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : basic_multiple_core_decomposition.h
* @brief      : a index method for temporal core decomposition
* @version    : 1.0
* @date       : 2021/08/09
******************************************************************************************************************/

#pragma once
#include "multiple_core_pair_map_index.h"
#include "multiple_core_pair_set_index.h"
#include "multiple_core_number_set_index.h"

namespace scnu{
    class branch_multiple_core_decomposition {
    public:

        static  void init(const shared_ptr<temporal_graph>& G,
                          const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                          const shared_ptr<thread_pool>& pool);

        static  void init(const shared_ptr<temporal_graph>& G,
                          const shared_ptr<unordered_map<uint32_t,
                                  shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                          const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                          const shared_ptr<thread_pool>& pool);

        static void init(const shared_ptr<temporal_graph> &G,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                         const shared_ptr<thread_pool>& pool);

        static  void init(const shared_ptr<temporal_graph>& G,
                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                          const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                          const shared_ptr<thread_pool>& pool);

        static uint32_t decompose(const shared_ptr<temporal_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                  const shared_ptr<thread_pool>& pool);

        static uint32_t decompose(const shared_ptr<temporal_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_index_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                  const shared_ptr<thread_pool>& pool);

        static uint32_t decompose(const shared_ptr<temporal_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                  const shared_ptr<thread_pool>& pool);

        static uint32_t decompose(const shared_ptr<temporal_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_index_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                  const shared_ptr<thread_pool>& pool);

    private:
        static void assign(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                           uint32_t k,
                           uint32_t h,
                           const shared_ptr<thread_pool>& pool);

        static shared_ptr<unordered_set<uint32_t>> find_core(const shared_ptr<temporal_graph>&G,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                             uint32_t k,
                                                             uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> find_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                             uint32_t k,
                                                             uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> find_core(const shared_ptr<temporal_graph>& G,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                             uint32_t k,
                                                             uint32_t h,
                                                             const shared_ptr<thread_pool>& pool);

        static shared_ptr<unordered_set<uint32_t>> find_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                             uint32_t k,
                                                             uint32_t h,
                                                             const shared_ptr<thread_pool>& pool);

        static shared_ptr<unordered_set<uint32_t>> find_h_core(const shared_ptr<temporal_graph>&G,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                               uint32_t k,
                                                               uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> find_h_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                               uint32_t k,
                                                               uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> find_h_core(const shared_ptr<temporal_graph>& G,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                               uint32_t k,
                                                               uint32_t h,
                                                               const shared_ptr<thread_pool>& pool);

        static shared_ptr<unordered_set<uint32_t>> find_h_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                               uint32_t k,
                                                               uint32_t h,
                                                               const shared_ptr<thread_pool>& pool);

        static shared_ptr<unordered_set<uint32_t>> find_k_core(const shared_ptr<temporal_graph>& G,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                               uint32_t k,
                                                               uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> find_k_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                               uint32_t k,
                                                               uint32_t h);

        static shared_ptr<unordered_set<uint32_t>> find_k_core(const shared_ptr<temporal_graph>& G,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                               uint32_t k,
                                                               uint32_t h,
                                                               const shared_ptr<thread_pool>& pool);

        static shared_ptr<unordered_set<uint32_t>> find_k_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                               uint32_t k,
                                                               uint32_t h,
                                                               const shared_ptr<thread_pool>& pool);

    };
}



