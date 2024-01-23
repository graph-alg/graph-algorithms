/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : wing_decomposition.h
* @brief      : Algorithms for k-wing decomposition
* @version    : 1.0
* @date       : 2020/09/06
******************************************************************************************************************/
#pragma once

#include "wing/wing_utility.h"
#include "wing/basic_wing_decomposition.h"

namespace scnu {
    /**
     * @brief a decomposition method
     * @cite Sarıyüce A E, Pinar A. Peeling bipartite networks for dense subgraph discovery[C]
     * //Proceedings of the Eleventh ACM International Conference on Web Search and Data Mining.
     * 2018: 504-512.
     */
    class peel_wing_decomposition {
    public:
        static uint32_t decompose(const std::shared_ptr<abstract_bipartite_graph> &G,
                                  const std::shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                  const std::shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>> &wing_edge_map);

    private:
        static shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>
        edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map);

        static shared_ptr<unordered_map<uint32_t, uint32_t>>
        vertex_priority_computation(const shared_ptr<abstract_bipartite_graph> &G);

    };
}


