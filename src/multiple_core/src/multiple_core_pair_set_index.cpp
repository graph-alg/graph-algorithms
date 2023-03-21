
#include "multiple_core/multiple_core_pair_set_index.h"

namespace scnu{
    multiple_core_pair_set_index::multiple_core_pair_set_index():
            core_pair_set(make_shared<unordered_set<pair<uint32_t,uint32_t>, hash_pair>>()){

    }

    bool multiple_core_pair_set_index::compare(const shared_ptr<multiple_core_pair_set_index>& other_multiple_core_pair_index){
        auto other_core_pair_set = other_multiple_core_pair_index->get_core_pair_set();
        return container_compare::same_associative_set(core_pair_set, other_core_pair_set);
    }

    bool multiple_core_pair_set_index::count(uint32_t k, uint32_t h) {
        return core_pair_set->count({k, h});
    }

    shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>> multiple_core_pair_set_index::get_core_pair_set(){
        return core_pair_set;
    }


    uint32_t multiple_core_pair_set_index::
    get_memory_cost()
    {
        uint32_t total_size = sizeof (*core_pair_set->begin()) * core_pair_set->size();
        return total_size;
    }


    void multiple_core_pair_set_index::insert(uint32_t k, uint32_t h){
        core_pair_set->insert({k, h});
    }
}

