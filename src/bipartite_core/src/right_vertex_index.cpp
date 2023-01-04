
#include "bipartite_core/right_vertex_index.h"

namespace scnu{
    right_vertex_index::right_vertex_index(): index_map(make_shared<map<uint32_t,uint32_t>>()) {

    }

    right_vertex_index::right_vertex_index(const shared_ptr<right_vertex_index>& other_right_index): right_vertex_index(){
        auto other_index_map = other_right_index->get_index_map();
        copy(other_index_map->begin(),other_index_map->end(),inserter(*index_map,index_map->end()));
    }

    bool right_vertex_index::compare(const shared_ptr<right_vertex_index>& another_right_node_index) {
        auto another_index_map = another_right_node_index->index_map;
        if (this->index_map->size() != another_index_map->size()) {
            return false;
        }

        for(const auto &[j,i]:*index_map){
            if(!another_index_map->count(j)){
                return false;
            }
            if(another_index_map->at(j)!=i){
                printf("%u,%u,%u\n",i,j,another_index_map->at(j));
                return false;
            }
        }
        return true;
    }

    void right_vertex_index::clear(){
        index_map->clear();
    }

    bool right_vertex_index::count(uint32_t j,uint32_t i) {
        return index_map->count(j) && index_map->at(j) >= i;
    }

    bool right_vertex_index::empty()
    {
        return index_map->empty();
    }

    shared_ptr<map<uint32_t , uint32_t>> right_vertex_index::get_index_map()
    {
        return index_map;
    }

    uint32_t right_vertex_index::get_i(uint32_t j)
    {
        if(!index_map->count(j)){
            return 0;
        }
        return index_map->at(j);
    }

    uint32_t right_vertex_index::get_maximal_j(uint32_t i)
    {
        for(auto iter = index_map->rbegin(); iter!=index_map->rend();++iter){
            if(iter->second >= i)
            {
                return iter->first;
            }
        }
        return 0;
    }

    void right_vertex_index::insert(uint32_t j, uint32_t i) {
        if (!index_map->count(j)) {
            index_map->insert({j, i});
        }else if (index_map->at(j) < i) {
            index_map->at(j) = i;
        }
    }

    void right_vertex_index::remove(uint32_t j){
        index_map->erase(j);
    }

    void right_vertex_index::remove(uint32_t j, uint32_t i) {
        if (!index_map->count(j)) {
            index_map->insert({j, i});
        }else if (index_map->at(j) > i) {
            index_map->at(j) = i;
        }
    }

    void right_vertex_index::set(uint32_t j, uint32_t i) {
        if (index_map->count(j)) {
            index_map->at(j) = i;
        }
    }

    bool right_vertex_index::operator==(const shared_ptr<right_vertex_index>& another_right_node_index) {
        return compare(another_right_node_index);
    }

    bool right_vertex_index::operator!=(const shared_ptr<right_vertex_index>& another_right_node_index) {
        return compare(another_right_node_index);
    }
}