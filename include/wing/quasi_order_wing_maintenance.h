/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : quasi_order_wing_maintenance.h
* @version    : 1.0
* @date       : 2021/05/11
******************************************************************************************************************/

#pragma once
#include "wing_utility.h"

namespace scnu{
    class quasi_order_wing_maintenance {
    public:

        static void insert(const shared_ptr<abstract_bipartite_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<
                                   extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_wing_order_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                           const shared_ptr<uint32_t>&
                           previous_k_max);


        static void insert(const shared_ptr<abstract_bipartite_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<
                                   extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_wing_order_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                           const shared_ptr<uint32_t>&previous_k_max,
                           uint32_t thread_count);

        static void remove(const shared_ptr<abstract_bipartite_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<
                                   extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                           const shared_ptr<uint32_t>& k_max);

        void remove(const shared_ptr<abstract_bipartite_graph> &G,
                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                    const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                    const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                    const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                    const shared_ptr<uint32_t> &k_max,
                    uint32_t thread_count);



    private:

        static void candidate_graph_finding(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_previous_k,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_k,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                            uint32_t k);

        static void candidate_graph_finding(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>> &edge_mutex_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_previous_k,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_k,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                            const shared_ptr<mutex> &F_k_mutex,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                            const shared_ptr<mutex> &evicted_edge_set_mutex,
                                            uint32_t k,
                                            uint32_t thread_count);

        static void find_affected_edge_set(const shared_ptr<abstract_bipartite_graph> &G,
                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &affected_edge_set,
                                           uint32_t k);

        static  void find_affected_edge_set(const shared_ptr<abstract_bipartite_graph> &G,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& affected_edge_set,
                                            const shared_ptr<mutex>& affected_edge_set_mutex,
                                            uint32_t k,
                                            uint32_t thread_count);

        static void init_insertion(const shared_ptr<abstract_bipartite_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS);

        static void init_insertion(const shared_ptr<abstract_bipartite_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map);

        static void partial_wing(const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& F_k,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                 uint32_t k);

        static uint32_t partial_wing_decomposition(const shared_ptr<abstract_bipartite_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                   uint32_t k);

        static uint32_t partial_wing_decomposition(const shared_ptr<abstract_bipartite_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list <uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                   uint32_t k,
                                                   uint32_t thread_count);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph>& G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& WS_k,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& F_k,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& evicted_edge_set,
                                             uint32_t k);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &WS_k,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &F_k,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                             const shared_ptr<mutex> &evicted_edge_set_mutex,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             uint32_t k,
                                             uint32_t thread_count);


        static void update_single_wing(const shared_ptr<abstract_bipartite_graph> &G,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                       uint32_t k);

        static void update_single_wing(const shared_ptr<abstract_bipartite_graph> &G,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                       const shared_ptr<mutex> &edge_set_mutex,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>>> &WS,
                                       uint32_t k,
                                       uint32_t thread_count);
    };
}



