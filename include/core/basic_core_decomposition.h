/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : basic_bipartite_core_decomposition.h
* @brief      : A basic core decomposition algorithm
* @version    : 1.0
* @date       : 2020/9/4
******************************************************************************************************************/
#pragma once
#include "core/core_utility.h"

namespace scnu
{
    /**
     * @details a class of basic core decomposition algorithms
     */
    class basic_core_decomposition
    {
    public:

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map);

        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  uint32_t thread_number);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<unordered_map<uint32_t,uint32_t>>>>& CD,
                                  uint32_t thread_number);


        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& core_order);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& core_order,
                                  uint32_t thread_number);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& tree);

        static uint32_t decompose(const shared_ptr<abstract_graph>& G,
                                  const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>>& tree,
                                  uint32_t thread_number);

        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                  const shared_ptr<gadget::Treap> &tree,
                                  const shared_ptr<vector<long long>> &root);

        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<double_list<uint32_t>>>>& k_order,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<double_node<uint32_t>>>> &node_map,
                                  const shared_ptr<gadget::Treap> &tree,
                                  const shared_ptr<vector<long long>> &root,
                                  uint32_t thread_number);

        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD);


        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                  uint32_t thread_number);

        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD);

        static uint32_t decompose(const shared_ptr<abstract_graph> &G,
                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,uint32_t>>>> &core_order,
                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_node<double,uint32_t>>>>& node_map,
                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>> &CD,
                                  uint32_t thread_number);

        static uint32_t get_core_degree(const shared_ptr<abstract_graph>& G,
                                        const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map,
                                        uint32_t v, uint32_t k);

        static void merger_set(const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &input_set_vector,
                               const shared_ptr<unordered_set<uint32_t>> &output_set,
                               const shared_ptr<thread_pool> &pool);
    };
}

