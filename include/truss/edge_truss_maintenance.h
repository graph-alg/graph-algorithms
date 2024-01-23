/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : edge_truss_maintenance.h
* @brief      : a k-truss maintenance method
* @version    : 1.0
* @date       : 2022/6/18
******************************************************************************************************************/
#pragma once
#include "truss_utility.h"

namespace scnu{
    class edge_truss_maintenance {
    public:
        static shared_ptr<unordered_set<uint32_t>> get_common_neighbor_set(const shared_ptr<abstract_graph>&G,
                                                                    uint32_t u,
                                                                    uint32_t v);

        void insert(const shared_ptr<abstract_graph>& G,
                    const shared_ptr<abstract_edge> &e,
                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map);

        void remove(const shared_ptr<abstract_graph>& G,
                    const shared_ptr<abstract_edge> &e,
                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map);

        pair<uint32_t, uint32_t> get_k1_k2(const shared_ptr<abstract_graph>& G,
                                           const shared_ptr<abstract_edge>&e1,
                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& edge_truss_map);

    };
}



