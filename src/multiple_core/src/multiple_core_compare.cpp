
#include "multiple_core/multiple_core_compare.h"


namespace scnu{
    bool multiple_core_compare::same(
            const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>> &core_map1,
            const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>> &core_map2) {
        return sub(core_map1, core_map2) && sub(core_map2, core_map1);
    }

    bool multiple_core_compare::same(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map1,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map2) {

        return sub(vertex_index_map1,vertex_index_map2) && sub(vertex_index_map2, vertex_index_map1);

    }

    bool multiple_core_compare::same(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_set_index>>> &vertex_index_map1,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_set_index>>> &vertex_index_map2) {

        return sub(vertex_index_map1,vertex_index_map2) && sub(vertex_index_map2, vertex_index_map1);
    }

    bool multiple_core_compare::same(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_number_set_index>>> &vertex_index_map1,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_number_set_index>>> &vertex_index_map2) {

        return sub(vertex_index_map1,vertex_index_map2) && sub(vertex_index_map2, vertex_index_map1);
    }

    bool multiple_core_compare::sub(
            const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>> &core_map1,
            const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<uint32_t>>, hash_pair>> &core_map2) {
        if(core_map1->size() !=core_map2->size()){
            printf("Unequal Size!\n");
            return false;
        }
        for(const auto&[p1, p1_set]:*core_map1){
            if(!core_map2->count(p1)){
                printf("Not Exist Pair: %u, %u\n", p1.first, p1.second);
                return false;
            }

            return container_compare::same_associative_set(p1_set, core_map2->at(p1));
        }
        return true;
    }

    bool multiple_core_compare::sub(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map1,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map2) {
        if(vertex_index_map1->size()!=vertex_index_map2->size())
        {
            printf("Unequal Size!\n");
            return false;
        }

        for(const auto&[v1,v1_index]:*vertex_index_map1){
            if(!vertex_index_map2->count(v1)){
                printf("Not Exist Vertex: %u\n", v1);
                return false;
            }
            auto v2_index = vertex_index_map2->at(v1);
            if(!v1_index->compare(v2_index)){
                printf("Unequal Vertex: %u\n", v1);
                return false;
            }
        }
        return true;
    }

    bool multiple_core_compare::sub(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_set_index>>> &vertex_index_map1,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_set_index>>> &vertex_index_map2) {
        if(vertex_index_map1->size()!=vertex_index_map2->size())
        {
            printf("Unequal Size!\n");
            return false;
        }

        for(const auto&[v1,v1_index]:*vertex_index_map1){
            if(!vertex_index_map2->count(v1)){
                printf("Not Exist Vertex!");
                return false;
            }
            auto v2_index = vertex_index_map2->at(v1);
            if(!v1_index->compare(v2_index)){
                printf("Unequal Vertex!");
                return false;
            }
        }
        return true;
    }

    bool multiple_core_compare::sub(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_number_set_index>>> &vertex_index_map1,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_number_set_index>>> &vertex_index_map2) {
        if(vertex_index_map1->size()!=vertex_index_map2->size())
        {
            printf("Unequal Size!\n");
            return false;
        }

        for(const auto&[v1,v1_index]:*vertex_index_map1){
            if(!vertex_index_map2->count(v1)){
                printf("Not Exist Vertex!");
                return false;
            }
            auto v2_index = vertex_index_map2->at(v1);
            if(!v1_index->compare(v2_index)){
                printf("Unequal Vertex!");
                return false;
            }
        }
        return true;
    }
}