/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : multiple_core_compare.h
* @brief      : a header file for temporal core compare
* @version    : 1.0
* @date       : 2021/09/20
******************************************************************************************************************/

#pragma once
#include "multiple_core/multiple_core_pair_map_index.h"
#include "multiple_core/multiple_core_pair_set_index.h"
#include "multiple_core/multiple_core_number_set_index.h"

namespace scnu{
    class multiple_core_compare {
    public:
        static bool same(const shared_ptr<unordered_map<pair<uint32_t, uint32_t >, shared_ptr<unordered_set<uint32_t>>, hash_pair>>& core_map1,
                         const shared_ptr<unordered_map<pair<uint32_t, uint32_t >,shared_ptr<unordered_set<uint32_t>>, hash_pair>>& core_map2);

        static bool same(const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map1,
                         const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map2);

        static bool same(const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_set_index>>>& vertex_index_map1,
                         const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_set_index>>>& vertex_index_map2);

        static bool same(const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_number_set_index>>>& vertex_index_map1,
                         const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_number_set_index>>>& vertex_index_map2);
    private:
        static bool sub(const shared_ptr<unordered_map<pair<uint32_t, uint32_t >, shared_ptr<unordered_set<uint32_t>>, hash_pair>>& core_map1,
                     const shared_ptr<unordered_map<pair<uint32_t, uint32_t >,shared_ptr<unordered_set<uint32_t>>, hash_pair>>& core_map2);

        static bool sub(const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map1,
                        const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map2);

        static bool sub(const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_set_index>>>& vertex_index_map1,
                        const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_set_index>>>& vertex_index_map2);

        static bool sub(const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_number_set_index>>>& vertex_index_map1,
                        const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_number_set_index>>>& vertex_index_map2);
    };
}





