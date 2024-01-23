/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : order_truss_decomposition.h
* @version    : 1.0
* @date       : 2022/06/18
******************************************************************************************************************/

#pragma once
#include "truss_utility.h"

namespace scnu{
    /**
     * @brief a order-based wing maintenance methods
     * @cite Zhang Y, Yu J X. Unboundedness and efficiency of truss maintenance in evolving graphs[C]
     * //Proceedings of the 2019 International Conference on Management of Data. 2019: 1024-1041.
     */
    class order_truss_maintenance {
    public:
        static void init(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext);

        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ts,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext);

        static void remove(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ts,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem);
    private:


        static uint32_t compute_edge_support(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<abstract_edge> &e_1,
                                             const shared_ptr<abstract_edge> &e_pivot,
                                             uint32_t k);

        static void edge_support_computation(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &s,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                             uint32_t k);

        static uint32_t k_new_computation(const shared_ptr<abstract_graph> &G,
                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                          const shared_ptr<abstract_edge> &e);

        static void remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &s,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                             const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &B,
                                             const shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>& e_pivot_node,
                                             const shared_ptr<abstract_edge>& e_current,
                                             uint32_t k);


        static bool
        test_order(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                   const shared_ptr<unordered_map<uint32_t ,shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                   const shared_ptr<abstract_edge> &e1,
                   const shared_ptr<abstract_edge> &e2);

        static bool
        test_order(const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &order_list,
                   const shared_ptr<abstract_edge> &e1,
                   const shared_ptr<abstract_edge> &e2);

        static void update_edge_support(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                        const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &B,
                                        const shared_ptr<abstract_edge> &e,
                                        uint32_t k);

        static void update_edge_support(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                        const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &B,
                                        uint32_t k);


        static void decrease_rem(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                 const shared_ptr<abstract_edge> &e,
                                 uint32_t k_new);

        static void increase_rem(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                const shared_ptr<abstract_edge> &e);

        static void update_rem_and_ts(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ts,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Q);


        static void update_order_list(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &s,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                      const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                      const shared_ptr<vector<shared_ptr<abstract_edge>>> &e_vector,
                                      uint32_t k);

        static void increase_ts(const shared_ptr<abstract_graph> &G,
                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ts,
                                uint32_t k);
    };
}






