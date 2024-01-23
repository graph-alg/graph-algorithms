/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bin_sort_truss_decomposition.h
* @brief      : a k-truss decomposition using bin sort
* @version    : 1.0
* @date       : 2022/6/18
******************************************************************************************************************/

#pragma once
#include "truss/truss_utility.h"

namespace scnu{
    class bin_sort_truss_decomposition {
    public:
        static void init(const shared_ptr<abstract_graph>& G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map);

        static void init(const shared_ptr<abstract_graph> &G,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map);

        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem);

        static void edge_support_computation(const shared_ptr<abstract_graph> &G,
                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_support_map);

    };
}



