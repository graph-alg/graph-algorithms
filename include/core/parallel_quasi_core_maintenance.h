/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : parallel_quasi_core_maintenance.h
* @brief      : Parallel core maintenance from graph view
* @version    : 1.0
* @date       : 2020/10/5
******************************************************************************************************************/

#pragma once
#include "core/core_utility.h"


namespace scnu {
    /**
     * @details core maintenance algorithms based on subgraph operations
     */
    class parallel_quasi_core_maintenance {
    public:

        static void
        bottom_up_insertion_maintenance(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                        uint32_t k_max,
                                        uint32_t thread_number);

        static void bottom_up_removal_maintenance(const shared_ptr<abstract_graph> &G,
                                                  const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &removal_edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                  uint32_t k_max,
                                                  uint32_t thread_number);

        static void initialize(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map);

        static void top_down_insertion_maintenance(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &insertion_edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &top_vertex_core_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                   uint32_t k_max,
                                                   long long n,
                                                   uint32_t thread_number);

        static void top_down_removal_maintenance(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &removal_edge_set,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &top_vertex_core_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                 uint32_t k_max,
                                                 long long n,
                                                 uint32_t thread_number);

    private:

        static shared_ptr<unordered_set<uint32_t>>
        bottom_up_candidate_graph(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                  uint32_t k,
                                  uint32_t thread_number);

        static void compute_estimated_core_number(const shared_ptr<abstract_graph> &G,
                                                  const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &ECD,
                                                  uint32_t thread_number);

        static uint32_t compute_estimated_core_number(const shared_ptr<abstract_graph> &G,
                                                      uint32_t w);


        static uint32_t edge_partition(const shared_ptr<abstract_graph> &G,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &core_edge_map);


        static shared_ptr<unordered_set<uint32_t>>
        top_down_candidate_graph(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &quasi_core_edge_map,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &top_vertex_core_map,
                                 uint32_t k,
                                 uint32_t thread_number);

        static void partial_core(const shared_ptr<unordered_set<uint32_t>> &Vc,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                 uint32_t k);

        static uint32_t partial_core_decomposition(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                   uint32_t k,
                                                   uint32_t thread_number);

        static shared_ptr<unordered_set<uint32_t>> remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &CDk,
                                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                                               const shared_ptr<unordered_set<uint32_t>> &evicted,
                                                                               uint32_t k,
                                                                               uint32_t thread_number);

        static void remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &CDk,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                const shared_ptr<unordered_set<uint32_t>> &evicted,
                                                uint32_t k,
                                                uint32_t thread_number);


        static shared_ptr<unordered_set<uint32_t>> update_single_core(const shared_ptr<abstract_graph> &G,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                                      uint32_t k,
                                                                      uint32_t thread_number);

        static shared_ptr<unordered_set<uint32_t>> update_single_core(const shared_ptr<abstract_graph> &G,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &quasi_core_edge_map,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                                      uint32_t k,
                                                                      uint32_t thread_number);

        static uint32_t vertex_quasi_core_decomposition(const shared_ptr<abstract_graph> &subgraph,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &degree_map,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &key_mutex_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_quasi_core_map,
                                                        const shared_ptr<thread_pool> &pool);
    };
}





