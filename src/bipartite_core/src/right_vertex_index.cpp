
#include "bipartite_core/bipartite_core_right_store_index.h"

namespace scnu{
    bipartite_core_right_store_index::bipartite_core_right_store_index(): index_map(make_shared<map<uint32_t,uint32_t>>()) {

    }

    bipartite_core_right_store_index::bipartite_core_right_store_index(const shared_ptr<bipartite_core_right_store_index>& other_right_index): bipartite_core_right_store_index(){
        auto other_index_map = other_right_index->get_index_map();
        copy(other_index_map->begin(),other_index_map->end(),inserter(*index_map,index_map->end()));
    }

    bool bipartite_core_right_store_index::compare(const shared_ptr<bipartite_core_right_store_index>& another_right_node_index) {
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

    void bipartite_core_right_store_index::clear(){
        index_map->clear();
    }

    bool bipartite_core_right_store_index::count(uint32_t j, uint32_t i) {
        return index_map->count(j) && index_map->at(j) >= i;
    }

    bool bipartite_core_right_store_index::empty()
    {
        return index_map->empty();
    }

    shared_ptr<map<uint32_t , uint32_t>> bipartite_core_right_store_index::get_index_map()
    {
        return index_map;
    }

    uint32_t bipartite_core_right_store_index::get_i(uint32_t j)
    {
        if(!index_map->count(j)){
            return 0;
        }
        return index_map->at(j);
    }

    uint32_t bipartite_core_right_store_index::get_maximal_j(uint32_t i)
    {
        for(auto iter = index_map->rbegin(); iter!=index_map->rend();++iter){
            if(iter->second >= i)
            {
                return iter->first;
            }
        }
        return 0;
    }

    void bipartite_core_right_store_index::insert(uint32_t j, uint32_t i) {
        if (!index_map->count(j)) {
            index_map->insert({j, i});
        }else if (index_map->at(j) < i) {
            index_map->at(j) = i;
        }
    }

    void bipartite_core_right_store_index::remove(uint32_t j){
        index_map->erase(j);
    }

    void bipartite_core_right_store_index::remove(uint32_t j, uint32_t i) {
        if (!index_map->count(j)) {
            index_map->insert({j, i});
        }else if (index_map->at(j) > i) {
            index_map->at(j) = i;
        }
    }

    void bipartite_core_right_store_index::set(uint32_t j, uint32_t i) {
        if (index_map->count(j)) {
            index_map->at(j) = i;
        }
    }

    bool bipartite_core_right_store_index::operator==(const shared_ptr<bipartite_core_right_store_index>& another_right_node_index) {
        return compare(another_right_node_index);
    }

    bool bipartite_core_right_store_index::operator!=(const shared_ptr<bipartite_core_right_store_index>& another_right_node_index) {
        return compare(another_right_node_index);
    }
}