
#include "bipartite_core/bipartite_core_left_store_index.h"

namespace scnu{
    bipartite_core_left_store_index::bipartite_core_left_store_index(): index_map(make_shared<map<uint32_t,uint32_t>>()) {

    }

    bipartite_core_left_store_index::bipartite_core_left_store_index(const shared_ptr<bipartite_core_left_store_index>& other_left_node_index): bipartite_core_left_store_index(){
        auto other_index_map = other_left_node_index->get_index_map();
        copy(other_index_map->begin(),other_index_map->end(),inserter(*index_map,index_map->end()));
    }

    bool bipartite_core_left_store_index::compare(shared_ptr<bipartite_core_left_store_index> &another_left_node_index) {
        auto another_index_map = another_left_node_index->index_map;
        if (index_map->size() != another_index_map->size()) {
            printf("Unequal Left Map Size!\n");
            return false;
        }

        for(const auto&[i,j]:*index_map){
            if(!another_index_map->count(i)){
                printf("Unequal Left Key!\n");
                return false;
            }
            if(another_index_map->at(i)!=j)
            {
                printf("%u,%u,%u",i,j,another_index_map->at(i));
                return false;
            }
        }
        return true;
    }

    void bipartite_core_left_store_index::clear(){
        index_map->clear();
    }

    bool bipartite_core_left_store_index::count(uint32_t i) {
        return index_map->count(i);
    }

    bool bipartite_core_left_store_index::count(uint32_t i, uint32_t j) {
        return index_map->count(i) && index_map->at(i) >= j;
    }

    bool bipartite_core_left_store_index::empty()
    {
        return index_map->empty();
    }

    shared_ptr<map<uint32_t , uint32_t>> bipartite_core_left_store_index::get_index_map()
    {
        return index_map;
    }

    uint32_t bipartite_core_left_store_index::get_j(uint32_t i)
    {
        if(!index_map->count(i)){
            return 0;
        }
        return index_map->at(i);
    }

    uint32_t bipartite_core_left_store_index::get_maximal_i(uint32_t j){
        for(auto iter = index_map->rbegin(); iter!=index_map->rend();++iter){
            if(iter->second >= j){
                return iter->first;
            }
        }
        return 0;
    }

    void bipartite_core_left_store_index::insert(uint32_t i, uint32_t j) {
        if (!index_map->count(i)) {
            index_map->insert({i, j});
        }else if (index_map->at(i) < j)
        {
            index_map->at(i) = j;
        }
    }

    void bipartite_core_left_store_index::remove(uint32_t i){
        index_map->erase(i);
    }

    void bipartite_core_left_store_index::remove(uint32_t i, uint32_t j) {
        if(!index_map->count(i)){
            index_map->insert({i, j});
        }
        else if(index_map->at(i) > j)
        {
            index_map->at(i) = j;
        }
    }

    void bipartite_core_left_store_index::set(uint32_t i, uint32_t j) {
        if(index_map->count(i)){
            index_map->at(i) = j;
        }
    }

    bool bipartite_core_left_store_index::operator==(shared_ptr<bipartite_core_left_store_index> &another_left_node_index) {
        return compare(another_left_node_index);
    }

    bool bipartite_core_left_store_index::operator!=(shared_ptr<bipartite_core_left_store_index> &another_left_node_index) {
        return !compare(another_left_node_index);
    }
}