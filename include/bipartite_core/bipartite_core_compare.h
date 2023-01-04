/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_core_compare.h
* @details    : a head file for bipartite core index compare
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core.h"
#include "bipartite_core/left_vertex_index.h"
#include "bipartite_core/right_vertex_index.h"

namespace scnu{
    class bipartite_core_compare {
    public:
        static void convert(
                const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map);

        static void convert(const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                            const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map);

        static void convert(const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                            const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                            const shared_ptr<unordered_set<pair<uint32_t, uint32_t>,hash_pair, equal_pair>> &core_set);

        static bool
        same(const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map1,
             const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map2);

        static bool same(const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map1,
             const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map1,
             const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map1);

        static bool
        sub(const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map1,
            const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map2);

        static bool same(const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map1,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map1,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map2,
                         const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map2);

        static bool sub(const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map1,
                        const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map1,
                        const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map2,
                        const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map2);
    };
}



