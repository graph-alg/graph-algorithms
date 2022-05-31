/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : order_maintenance.h
* @brief      : An core-order based maintenance algorithm
* @version    : 1.1
* @date       : 2020/10/20
******************************************************************************************************************/

#pragma once

#include "core/core_utility.h"
#include "core/quasi_core_maintenance.h"

namespace scnu {
    /**
     * @details An core-order based maintenance algorithm, it can tackle a set of edges simultaneously
     * @remarks core-order is strictly ordered
     * @cite  Zhang Y, Yu J X. Unboundedness and efficiency of truss maintenance in evolving graphs[C]
     * Proceedings of the 2019 International Conference on Management of Data. 2019: 1024-1041.
     */
    class order_maintenance {
    public:
        static void batch_insert(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &rem,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &ext
        );

        static void batch_insert(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                 const shared_ptr<unordered_map<uint32_t ,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>> &CD,
                                 uint32_t previous_max_k);

        static void batch_remove(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ER,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &rem,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &ts
        );

        static void init_order(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<
                                       extend_list<double, uint32_t>>>> &k_order,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &rem,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &ext
        );

    private:
        static uint32_t get_core_degree(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                        uint32_t k,
                                        uint32_t u);

        static uint32_t recompute_core_number(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                              uint32_t u);

        static bool test_order(const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                               uint32_t u, uint32_t
                               v);

        static shared_ptr<extend_node<double, uint32_t>> partial_core(const shared_ptr<abstract_graph> &G,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                                                      const shared_ptr<unordered_set<uint32_t>> &Ck,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &s,
                                                                      uint32_t k);

        static shared_ptr<extend_node<double, uint32_t>> partial_core(const shared_ptr<abstract_graph> &G,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                                                      const shared_ptr<unordered_set<uint32_t>> &Ck,
                                                                      const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                                      uint32_t k);

        static void partial_core(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_set<uint32_t>> &Vc,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_node<double, uint32_t>>>> &node_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                 const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                 uint32_t k);

        static void partial_core_decomposition(const shared_ptr<abstract_graph>& G,
                                               const shared_ptr<unordered_set<uint32_t>>& Vc,
                                               const shared_ptr<unordered_set<uint32_t>>& evicted_set,
                                               const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                               const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD,
                                               const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                               const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                               uint32_t k);


        static void remove_unsatisfied_vertices(const shared_ptr<abstract_graph>& G,
                                                const shared_ptr<unordered_set<uint32_t>>& Vc,
                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD,
                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                                const shared_ptr<unordered_set<uint32_t>>& evicted_set,
                                                uint32_t w,
                                                uint32_t k);

        static void remove_unsatisfied_vertices(const shared_ptr<abstract_graph>& G,
                                                const shared_ptr<unordered_set<uint32_t>>& Ck,
                                                uint32_t w,
                                                uint32_t k,
                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t ,uint32_t>>>>& CD,
                                                const shared_ptr<unordered_map<uint32_t ,shared_ptr<extend_list<double,uint32_t>>>>& k_order,
                                                const shared_ptr<unordered_set<uint32_t>>& evicted_set);

        static void remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                const shared_ptr<unordered_set<uint32_t>> &evictedSet,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &kMcd,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &kOrder,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_node<double, uint32_t>>>> &nodeMap,
                                                uint32_t k);
        static void search_influenced_vertices(const shared_ptr<abstract_graph> &G,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_node<double, uint32_t>>>> &node_map,
                                               const shared_ptr<unordered_set<uint32_t>> &Vc,
                                               const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                               uint32_t k);

        static void search_influenced_vertices(const shared_ptr<abstract_graph> &G,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &k_order,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                               const shared_ptr<unordered_set<uint32_t>> &Ck,
                                               const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                               const shared_ptr<extend_node<double,uint32_t>> &old_head,
                                               uint32_t k);
    };
}



