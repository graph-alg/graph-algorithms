
#include "multiple_core/multiple_core_number_set_index.h"

namespace scnu{
    multiple_core_number_set_index::multiple_core_number_set_index() :
            core_number_map(make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair>>()) {

    }


    bool multiple_core_number_set_index::compare(const shared_ptr<multiple_core_number_set_index>& other_multiple_core_set_index){
        auto other_multiple_core_map = other_multiple_core_set_index->get_multiple_core_map();

        auto result = true;
        if(core_number_map->size() != other_multiple_core_map->size()){
            result =  false;
        }
        for(const auto&p:*core_number_map){
            if(!other_multiple_core_map->count(p))
            {
                result = false;
                break;
            }
        }
        return result;
    }

    bool multiple_core_number_set_index::count(uint32_t k, uint32_t h)
    {
        for(const auto &p:*core_number_map){
            if(p.first >= k && p.second >= h){
                return true;
            }
        }
        return false;
    }

    shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>>
    multiple_core_number_set_index::get_multiple_core_map() {
        return core_number_map;
    }

    uint32_t multiple_core_number_set_index::get_memory_cost()
    {
        return sizeof(*core_number_map->begin()) * core_number_map->size();
    }


    void multiple_core_number_set_index::insert(uint32_t k, uint32_t h) {
        if (k == 0 || h == 0) {
            return;
        }
        for(auto iter = core_number_map->begin(); iter!=core_number_map->end();){
            auto &p = *iter;
            ++iter;
            if((p.first <= k && p.second < h) || (p.first < k && p.second <= h)){
                core_number_map->erase(p);
            }
        }
        core_number_map->insert({k, h});
    }
}