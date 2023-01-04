
#include "bipartite_core/left_vertex_index.h"

namespace scnu{
    left_vertex_index::left_vertex_index(): index_map(make_shared<map<uint32_t,uint32_t>>()) {

    }

    left_vertex_index::left_vertex_index(const shared_ptr<left_vertex_index>& other_left_node_index): left_vertex_index(){
        auto other_index_map = other_left_node_index->get_index_map();
        copy(other_index_map->begin(),other_index_map->end(),inserter(*index_map,index_map->end()));
    }

    bool left_vertex_index::compare(shared_ptr<left_vertex_index> &another_left_node_index) {
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

    void left_vertex_index::clear(){
        index_map->clear();
    }

    bool left_vertex_index::count(uint32_t i) {
        return index_map->count(i);
    }

    bool left_vertex_index::count(uint32_t i, uint32_t j) {
        return index_map->count(i) && index_map->at(i) >= j;
    }

    bool left_vertex_index::empty()
    {
        return index_map->empty();
    }

    shared_ptr<map<uint32_t , uint32_t>> left_vertex_index::get_index_map()
    {
        return index_map;
    }

    uint32_t left_vertex_index::get_j(uint32_t i)
    {
        if(!index_map->count(i)){
            return 0;
        }
        return index_map->at(i);
    }

    uint32_t left_vertex_index::get_maximal_i(uint32_t j){
        for(auto iter = index_map->rbegin(); iter!=index_map->rend();++iter){
            if(iter->second >= j){
                return iter->first;
            }
        }
        return 0;
    }

    void left_vertex_index::insert(uint32_t i, uint32_t j) {
        if (!index_map->count(i)) {
            index_map->insert({i, j});
        }else if (index_map->at(i) < j)
        {
            index_map->at(i) = j;
        }
    }

    void left_vertex_index::remove(uint32_t i){
        index_map->erase(i);
    }

    void left_vertex_index::remove(uint32_t i, uint32_t j) {
        if(!index_map->count(i)){
            index_map->insert({i, j});
        }
        else if(index_map->at(i) > j)
        {
            index_map->at(i) = j;
        }
    }

    void left_vertex_index::set(uint32_t i, uint32_t j) {
        if(index_map->count(i)){
            index_map->at(i) = j;
        }
    }

    bool left_vertex_index::operator==(shared_ptr<left_vertex_index> &another_left_node_index) {
        return compare(another_left_node_index);
    }

    bool left_vertex_index::operator!=(shared_ptr<left_vertex_index> &another_left_node_index) {
        return !compare(another_left_node_index);
    }
}