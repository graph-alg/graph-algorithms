/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : wing_decomposition.h
* @brief      : A common head files for data structure
* @version    : 1.0
* @date       : 2020/10/29
******************************************************************************************************************/

#pragma once
#include "peel_wing_decomposition.h"

namespace scnu
{
    class order_wing_maintenance_back
    {
    public:

        static void insertion_maintenance(const std::shared_ptr<abstract_bipartite_graph>& current_graph,
                                                   const std::shared_ptr<abstract_bipartite_graph>& insertion_graph,
                                                   const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                                                   const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rem_map
        );

        static void insertion_parallel_maintenance(const std::shared_ptr<abstract_bipartite_graph> &current_graph,
                                                               const std::shared_ptr<abstract_bipartite_graph> &insertion_graph,
                                                               const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_truss_order_map,
                                                               const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rem_map,
                                                               const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rectangle_map,
                                                               uint32_t thread_number
        );

        static void insertion_thread_maintenance(const std::shared_ptr<abstract_bipartite_graph>& current_graph,
                                                 const std::shared_ptr<abstract_bipartite_graph>& insertion_graph,
                                                 const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                                                 const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rem_map,
                                                 const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_truss_map
        );

        static void removel_order_maintenance(const std::shared_ptr<abstract_bipartite_graph> &previous_graph,
                                        const std::shared_ptr<abstract_bipartite_graph> &current_graph,
                                        const std::shared_ptr<abstract_bipartite_graph> &removel_graph,
                                        const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                                        const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rem_map,
                                        const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& previous_edge_ts_map
                );

        static void removel_parallel_maintenance(const std::shared_ptr<abstract_bipartite_graph> &previous_graph,
                                                 const std::shared_ptr<abstract_bipartite_graph> &removal_graph,
                                                 const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_truss_order_map,
                                                 const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>>>& WS,
                                                 int k_max,
                                                 uint32_t thread_number
        );

        static void removel_maintenance(const std::shared_ptr<abstract_bipartite_graph> &previous_graph,
                                        const std::shared_ptr<abstract_bipartite_graph> &removel_graph,
                                        const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_truss_order_map,
                                        const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>>>& WS,
                                        int k_max
        );


        static int order_truss_decomposition(const std::shared_ptr<abstract_bipartite_graph>& graph,
                                             const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edgeTrussMap,
                                             const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& trussEdgeMap,
                                             const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                                             const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rem_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>>> &wing_edge_support_map
        );
        static void get_all_edge_ts_num(const std::shared_ptr<abstract_bipartite_graph>& graph,
                                        const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                                        const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& previous_edge_ts_map
        );

        static void get_all_edge_rectangle_num(const std::shared_ptr<abstract_bipartite_graph> &graph,
                                                                const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_truss_order_map,
                                                                const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &previous_edge_ts_map
        );


    private:
        static shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>
        get_all_order(
                const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                uint32_t k
        );
        static int removal_quasi_truss_decomposition(const std::shared_ptr<abstract_bipartite_graph>& sub_graph,
                                                     const std::shared_ptr<abstract_bipartite_graph>& entire_graph,
                                                     const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& quasi_edge_truss_map,
                                                     const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& quasi_truss_edge_map,
                                                     const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>>>& WS
        );

        static void remove_parallel_edge_set(const std::shared_ptr<abstract_bipartite_graph>& G,
                                             const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>&quasi_truss_edge_map,
                                             const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>>>& WS,
                                             const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_truss_order_map,
                                             const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edgeTrussMap,
                                             int k,
                                             uint32_t thread_number
        );

        static void remove_edge_set(const std::shared_ptr<abstract_bipartite_graph>& G,
                                    const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>&quasi_truss_edge_map,
                                    const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>>>& WS,
                                    const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t, shared_ptr<abstract_bipartite_edge>>>>> &previous_truss_order_map,
                                    const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edgeTrussMap,
                                    int k);

        static int get_edge_ts_num(const std::shared_ptr<abstract_bipartite_graph>& graph,
                                    const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                                    const std::shared_ptr<abstract_bipartite_edge>& edge,
                                    uint32_t k
        );

        static void removal_unsatisfied_edges(const std::shared_ptr<abstract_bipartite_graph>&graph,
                                              const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rectangle_map,
                                              const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& rectangle_edge_map,
                                              const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_truss_map,
                                              const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& truss_edge_map,
                                              const std::shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<uint32_t,shared_ptr<abstract_bipartite_edge>>>>>& previous_truss_order_map,
                                              const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rem_map,
                                              const shared_ptr<unordered_map<uint32_t,
                                                      shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>>> &wing_edge_support_map,
                                              uint32_t k);

        static void removal_unsatisfied_edges(const std::shared_ptr<abstract_bipartite_graph>&graph,
                                                               const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rectangle_map,
                                                               const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& rectangle_edge_map,
                                                               const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_truss_map,
                                                               const std::shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& truss_edge_map,
                                                               uint32_t k);

        static uint32_t get_edge_rectangle_count(const shared_ptr<abstract_bipartite_graph>& G,
                                                 const shared_ptr<abstract_bipartite_edge>& e1);

        static void get_edge_rectangle_map(const shared_ptr<abstract_bipartite_graph>& G,
                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_rectangle_map,
                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>& rectangle_edge_map);


    };

}
