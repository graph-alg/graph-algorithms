/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : multiple_core_query.h
* @brief      : query algorithms
* @version    : 1.0
* @date       : 2022/03/13
******************************************************************************************************************/

#pragma  once
#include "multiple_core_pair_map_index.h"
#include "multiple_core_pair_set_index.h"
#include "multiple_core_number_set_index.h"

namespace scnu{
    class multiple_core_query {
    public:
        template<class index_container_type>
        static void query_multiple_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<index_container_type>>> &vertex_index_map,
                                        const shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>>& core_pair_set,
                                        const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>>& core_pair_map)
        {

            for(const auto &[k,h]:*core_pair_set){
                core_pair_map->insert({{k,h}, make_shared<unordered_set<uint32_t>>()});
                for(const auto&[u,u_index]:*vertex_index_map){
                    if(u_index->count(k, h)){
                        core_pair_map->at({k,h})->insert(u);
                    }
                }
            }
        }
    };
}


