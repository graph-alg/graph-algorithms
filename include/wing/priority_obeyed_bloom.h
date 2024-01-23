/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : priority_obeyed_bloom.h
* @brief      : a data structure for storing bloom
* @version    : 1.0
* @date       : 2021/4/6
******************************************************************************************************************/

#pragma once
#include "wing_utility.h"

namespace scnu {
    class priority_obeyed_bloom {
    public:

        priority_obeyed_bloom(uint32_t u, uint32_t w);

        uint32_t count(const shared_ptr<abstract_bipartite_edge>& e);

        bool empty();

        uint32_t get_butterfly_count() const;

        shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<abstract_bipartite_edge>>>
        get_edge_map();

        uint32_t get_k() const;

        shared_ptr<abstract_bipartite_edge> get_twin(const shared_ptr<abstract_bipartite_edge> &e);

        pair<uint32_t,uint32_t> get_vertex_pair();

        void remove_edge(const shared_ptr<abstract_bipartite_edge>&e);

        void set_butterfly_count(uint32_t value);

        void link_twin(const shared_ptr<abstract_bipartite_edge> &e,
                       const shared_ptr<abstract_bipartite_edge> &twin_edge);

    private:
        pair<uint32_t,uint32_t> vertex_pair;
        /**
         * @brief the total butterflies of this bloom
         */
        uint32_t butterfly_count;

        /**
         * @brief a map of edge and its twin
         */
        shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<abstract_bipartite_edge>>> edge_twin_map;
    };
}
