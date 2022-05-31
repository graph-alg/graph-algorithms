/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : quasi_core_maintenance.h
* @brief      : algorithms for k-core maintenance
* @version    : 1.2
* @date       : 2022/03/31
******************************************************************************************************************/

#pragma once
#include "core/core_utility.h"

namespace scnu
{
    /**
     * @details a graph-based core maintenance method, it divides maintenance process to several steps: (1) quasi-k-core
     * (2) candidate graph (3) partial-k-core
     */
    class quasi_core_maintenance
    {
    public:
        static void
        insert(const shared_ptr<abstract_graph> &previous_graph,
               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &insertion_edge_set,
               const shared_ptr<unordered_map<uint32_t,uint32_t>> &vertex_core_map,
               const shared_ptr<uint32_t>& k_max);


        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                           const shared_ptr<uint32_t> &k_max,
                           uint32_t thread_number);

        static void remove(const shared_ptr<abstract_graph> &previous_graph,
                    const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &removal_edge_set,
                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                    const shared_ptr<uint32_t>& k_max);


        static void remove(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &removal_edge_set,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                           const shared_ptr<uint32_t> &k_max,
                           uint32_t thread_number);


    private:

        static void candidate_graph_finding(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &quasi_core_edge_map,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>>& CD,
                                            const shared_ptr<unordered_set<uint32_t>>& Vc,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>>& Vc_degree_map,
                                            uint32_t k);


        static void candidate_graph_finding(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> & CD,
                                            const shared_ptr<unordered_set<uint32_t>> &V,
                                            const shared_ptr<unordered_set<uint32_t>> &Vc,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &Vc_degree_map,
                                            const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& current_task_vector,
                                            const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& evicted_task_vector,
                                            uint32_t k,
                                            const shared_ptr<thread_pool> &pool);


        static  uint32_t get_core_degree(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                         uint32_t u,
                                         uint32_t k);

        static void merger_set(const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &input_set_vector,
                               const shared_ptr<unordered_set<uint32_t>> &output_set,
                               const shared_ptr<thread_pool> &pool);


        static uint32_t partial_core_decomposition(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &CD,
                                                   uint32_t k);

        static uint32_t partial_core_decomposition(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_vector,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                   const shared_ptr<unordered_set<uint32_t>> &F_k,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                   const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &current_task_vector,
                                                   uint32_t k,
                                                   const shared_ptr<thread_pool>& pool);

        static uint32_t quasi_core_edge_decomposition(const shared_ptr<abstract_graph> &entire_graph,
                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                      const shared_ptr<unordered_map<uint32_t,
                                                                                shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>
                                                       &quasi_core_edge_map);


        static void remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                uint32_t w,
                                                uint32_t k,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &Vc_degree_map,
                                                const shared_ptr<unordered_set<uint32_t>>& evicted_set);

        static void remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                uint32_t k,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &Vc_degree_map);


        static void remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                const shared_ptr<unordered_map<uint32_t,uint32_t>> &Vc_degree_map,
                                                const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& current_task_vector,
                                                uint32_t k,
                                                const shared_ptr<thread_pool>& pool);


        static void update_single_core(const shared_ptr<abstract_graph> &G,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &quasi_core_edge_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &CD,
                                       const shared_ptr<unordered_set<uint32_t>> &previous_removed_vertex_set,
                                       const shared_ptr<unordered_set<uint32_t>> &current_removed_vertex_set,
                                       uint32_t k);

        static void update_single_core(const shared_ptr<abstract_graph> &G,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_vector,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &current_vertex_core_map,
                                       const shared_ptr<unordered_set<uint32_t>> &current_removed_vertex_set,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &CD,
                                       const shared_ptr<unordered_set<uint32_t>> &V,
                                       uint32_t k,
                                       const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &current_task_vector,
                                       const shared_ptr<thread_pool> &pool);

        static void update_single_core(const shared_ptr<abstract_graph> &G,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_vector,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                       const shared_ptr<unordered_set<uint32_t>> &current_removed_vertex_set,
                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &CD,
                                       const shared_ptr<unordered_set<uint32_t>> &V,
                                       uint32_t k,
                                       const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &current_task_vector,
                                       const shared_ptr<thread_pool> &pool);

    };
}