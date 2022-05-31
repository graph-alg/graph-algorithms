/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : order_core_maintenance.h
* @brief      : Order-based core maintenance methods
* @version    : 1.1
* @date       : 2020/10/16
******************************************************************************************************************/

#pragma once

#include "core/core_utility.h"
#include "core/traversal_core_maintenance.h"

namespace scnu
{
    /**
     * @class order_core_maintenance
     * @details a maintenance algorithm utilizing the removing order among decomposition process
     * @cite A fast order-based approach for core maintenance
     */
    class order_core_maintenance
    {
    public:
        static void insert(const shared_ptr<abstract_graph>&G,
                           const shared_ptr<abstract_edge>& e,
                           const shared_ptr<unordered_map<uint32_t,uint32_t>>& core,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>>& k_order,
                           const shared_ptr<unordered_map<uint32_t , shared_ptr<double_node<uint32_t>>>>& node_map,
                           const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> & tree,
                           const shared_ptr<unordered_map<uint32_t , uint32_t>>& remaining_degree_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_degree_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t> >>>& rcd,
                           uint32_t n);

        static void insert(const shared_ptr<abstract_graph> &G,
                           const shared_ptr<abstract_edge> &e,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>> &k_order,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                           const shared_ptr<gadget::Treap> &tree,
                           const shared_ptr<vector<long long>> &root,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map,
                           const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t> >>> &rcd,
                           uint32_t n
        );

        static void order_remove(const shared_ptr<abstract_graph> &G,
                                 const shared_ptr<abstract_edge> &e,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>> &core_vertex_order,
                                 const shared_ptr<unordered_map<uint32_t , shared_ptr<double_node<uint32_t>>>> &node_map,
                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> & core_order_query_map,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& rcd,
                                 uint32_t n);

        static  void remove(const shared_ptr<abstract_graph> &G,
                            const shared_ptr<abstract_edge> &e,
                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<double_list<uint32_t>>>> &core_vertex_order,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                            const shared_ptr<gadget::Treap> &tree,
                            const shared_ptr<vector<long long>> &root,
                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                            uint32_t n);

        static void initialize(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, uint32_t>>>> &tree,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map,
                               const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& mcd,
                               uint32_t n);

        static void initialize(const shared_ptr<abstract_graph> &G,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                               const shared_ptr<gadget::Treap> &tree,
                               const shared_ptr<vector<long long>> &root,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map,
                               const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                               uint32_t n);

    private:

        static uint32_t compute_rcd(const shared_ptr<abstract_graph> &G,
                                                            uint32_t u,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                            uint32_t h);

        static void multi_hop_prepare_rcd_insertion(const shared_ptr<abstract_graph> &G,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                    uint32_t n,
                                                    const shared_ptr<abstract_edge> &e);

        static void multi_hop_prepare_rcd_removal(const shared_ptr<abstract_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                                  uint32_t n,
                                                  const shared_ptr<abstract_edge> &e);

        static void multi_hop_recompute_rcd(const shared_ptr<abstract_graph> &G,
                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                            uint32_t n,
                                            const shared_ptr<list<uint32_t>> &changed,
                                            bool state);

        static void propagate_dismissal(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &rcd,
                                        const shared_ptr<unordered_map<uint32_t, long long>> &cd,
                                        const shared_ptr<unordered_set<uint32_t>> &dismissed,
                                        const shared_ptr<unordered_set<uint32_t>> &visited,
                                        uint32_t k,
                                        uint32_t v,
                                        uint32_t n,
                                        const shared_ptr<list<uint32_t>> &changed);


        static void remove_candidates(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<list<uint32_t>> &Vc,
                                      const shared_ptr<unordered_set<uint32_t>> &Vc_set,
                                      const shared_ptr<unordered_map<uint32_t,shared_ptr<double_node<uint32_t>>>> &node_map,
                                      const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &tree,
                                      const shared_ptr<map<double, shared_ptr<double_node<uint32_t>>>> &B,
                                      const shared_ptr<double_list<uint32_t>> &new_k_list,
                                      const shared_ptr<vector<uint32_t>> &swap,
                                      uint32_t w,
                                      uint32_t K,
                                      const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                                      const shared_ptr<unordered_map<uint32_t,uint32_t>> &remaining_degree_map,
                                      const shared_ptr<unordered_map<uint32_t,uint32_t>> &candidate_degree_map
        );

        static void remove_candidates(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<list<uint32_t>> &Vc,
                                      const shared_ptr<unordered_set<uint32_t>> &Vc_set,
                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                      const shared_ptr<gadget::Treap> &tree,
                                      const shared_ptr<vector<long long>> &root,
                                      const shared_ptr<map<double, shared_ptr<double_node<uint32_t>>>> &B,
                                      const shared_ptr<double_list<uint32_t>> &new_k_list,
                                      const shared_ptr<vector<uint32_t>> &swap,
                                      uint32_t w,
                                      uint32_t K,
                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &remaining_degree_map,
                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_degree_map
        );


        static bool test_order(const shared_ptr<unordered_map<uint32_t,uint32_t>> &core,
                               const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &tree,
                               uint32_t u,
                               uint32_t v);

        static bool test_order(const shared_ptr<unordered_map<uint32_t, uint32_t>> &core,
                               const shared_ptr<gadget::Treap> &tree,
                               uint32_t u,
                               uint32_t v);
    };
}



