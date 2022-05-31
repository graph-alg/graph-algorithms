/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : jes_core_maintenance.h
* @brief      : Generalized traversal methods
* @version    : 1.0
* @date       : 2020/10/21
******************************************************************************************************************/

#pragma once
#include "core/core_utility.h"

namespace scnu
{
    /**
     * @details A generalized traversal methods, it supports multi-threads
     * @cite Hua Q S, Shi Y, Yu D, et al. Faster Parallel bipartite_core Maintenance Algorithms in Dynamic Graphs[J].
     * IEEE Transactions on Parallel and Distributed Systems, 2019, 31(6): 1287-1300.
     */
    class jes_core_maintenance
    {
    public:

        static void joint_delete(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                 uint32_t thread_number);

        static void joint_insert(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                 uint32_t thread_number);

    private:

        static uint32_t compute_core_number(const shared_ptr<abstract_graph>& G,
                                       const shared_ptr<unordered_map<uint32_t,uint32_t>>& core,
                                       uint32_t v);

        static pair<shared_ptr<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>,
                shared_ptr<unordered_set<uint32_t>>>
        compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                const shared_ptr<thread_pool>& pool);

        static  pair<shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>,
        shared_ptr<unordered_set<uint32_t>>>
        compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                const shared_ptr<thread_pool>& pool);

        static uint32_t compute_pre_core_number(const shared_ptr<abstract_graph>& G,
                                                const shared_ptr<unordered_map<uint32_t,uint32_t>>& core,
                                                uint32_t v);

        static shared_ptr<unordered_set<uint32_t>> 
        get_maximal_three_hop_independent_set( const shared_ptr<abstract_graph> &G,
                                               const shared_ptr<unordered_set<uint32_t>> &superior_vertex_set);

        static  shared_ptr<unordered_set<uint32_t>>
        get_k_path_tree(const shared_ptr<abstract_graph> &G,
                        const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                        uint32_t w,
                        uint32_t k);

        static shared_ptr<unordered_set<uint32_t>>
        get_k_path_tree(const shared_ptr<abstract_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                  const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                  uint32_t k
        );


        static uint32_t
        get_superior_degree(const shared_ptr<abstract_graph> &G,
                            const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                            uint32_t u);

        static shared_ptr<unordered_set<uint32_t>>
        get_superior_vertex_set(const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &E,
                                const shared_ptr<unordered_map<uint32_t,uint32_t>> &core);

        static shared_ptr<unordered_set<uint32_t>>
        k_joint_delete(const shared_ptr<abstract_graph>& G,
                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& Ek,
                       const shared_ptr<unordered_map<uint32_t,uint32_t>>& core,
                       uint32_t k);

        static shared_ptr<unordered_set<uint32_t>>
        k_joint_insert(const shared_ptr<abstract_graph> &G,
                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                       const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                       uint32_t k);
    };
}




