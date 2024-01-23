/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : index_wing_decomposition.h
* @brief      : an index decomposition
* @version    : 1.0
* @date       : 2021/4/6
******************************************************************************************************************/

#pragma once
#include "wing_utility.h"
#include "priority_obeyed_bloom.h"

namespace scnu{
    class BE_index{
    public:
        BE_index();

        bool count(uint32_t u,uint32_t v);

        bool count(const shared_ptr<abstract_bipartite_edge>& e);

        uint32_t get_support(const shared_ptr<abstract_bipartite_edge>& e);

        shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>> get_edge_support_map();

        shared_ptr<priority_obeyed_bloom> get_bloom(uint32_t u,uint32_t v);

        shared_ptr<unordered_map<pair<uint32_t,uint32_t>,shared_ptr<priority_obeyed_bloom>,hash_pair,equal_pair>>
        get_bloom_map();

        shared_ptr<unordered_set<shared_ptr<priority_obeyed_bloom>>>
        get_bloom_set(const shared_ptr<abstract_bipartite_edge>& e);

        void insert_bloom(const shared_ptr<priority_obeyed_bloom>& B);

        void link_bloom(const shared_ptr<abstract_bipartite_edge>& e,
                        const shared_ptr<priority_obeyed_bloom>& B);

        void insert_edge(const shared_ptr<abstract_bipartite_edge>& e,
                         uint32_t support);

        void removal_bloom(const shared_ptr<priority_obeyed_bloom>& B);

        void remove_edge(const shared_ptr<abstract_bipartite_edge>& e);

        void remove_edge(const shared_ptr<abstract_bipartite_edge>& e,
                         const shared_ptr<priority_obeyed_bloom>& B);

        void update_support(const shared_ptr<abstract_bipartite_edge>& e,
                                      uint32_t value);

    private:
        shared_ptr<unordered_map<pair<uint32_t,uint32_t>,shared_ptr<priority_obeyed_bloom>,hash_pair,equal_pair>> bloom_map;
        shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>> edge_support_map;
        shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<unordered_set<shared_ptr<priority_obeyed_bloom>>>>> edge_bloom_map;
    };
}
