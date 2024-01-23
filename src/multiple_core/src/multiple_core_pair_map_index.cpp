
#include "multiple_core/multiple_core_pair_map_index.h"

namespace scnu{
    multiple_core_pair_map_index::multiple_core_pair_map_index()
                : left_map(make_shared<unordered_map<uint32_t, uint32_t>>()), right_map(make_shared<unordered_map<uint32_t, uint32_t>>())  {
    }

    multiple_core_pair_map_index::multiple_core_pair_map_index(const shared_ptr<multiple_core_pair_map_index>& other_multiple_core_map_index)
                                                            : multiple_core_pair_map_index()
    {
        for(const auto &[h, k]:*other_multiple_core_map_index->left_map){
            left_map->insert({h, k});
        }
        for(const auto &[k,h]:*other_multiple_core_map_index->right_map){
            right_map->insert({k, h});
        }
    }

    void multiple_core_pair_map_index::clear(){
        left_map->clear();
        right_map->clear();
    }

    bool multiple_core_pair_map_index::compare(const shared_ptr<multiple_core_pair_map_index>& other_multiple_core_map_index){
        auto other_left_map = other_multiple_core_map_index->get_left_map();
        auto other_right_map = other_multiple_core_map_index->get_right_map();

        return container_compare::same_associative_map(left_map, other_left_map)
                && container_compare::same_associative_map(right_map, other_right_map);
    }

    bool multiple_core_pair_map_index::count(uint32_t k, uint32_t h)
    {
        return k > h ? (left_map->count(h) && left_map->at(h) >=  k) : (right_map->count(k) && right_map->at(k) >= h);
    }

    bool multiple_core_pair_map_index::empty(){
        return left_map->empty() && right_map->empty();
    }

    bool multiple_core_pair_map_index::equal(uint32_t k, uint32_t h) {
        if(k > h){
            return  left_map->count(h) && left_map->at(h) == k;
        }else if(k < h){
            return  right_map->count(k) && right_map->at(k) == h;
        }else{
            return  (left_map->count(h) && left_map->at(h) == k) || (right_map->count(k) && right_map->at(k) == h);
        }
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>> multiple_core_pair_map_index::get_left_map(){
        return left_map;
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>> multiple_core_pair_map_index::get_right_map(){
        return right_map;
    }

    uint32_t multiple_core_pair_map_index::get_memory_cost()
    {
        return sizeof(*left_map->begin()) * left_map->size()
               +sizeof(*right_map->begin()) * right_map->size();
    }

    uint32_t multiple_core_pair_map_index::get_delta(){
        return left_map->size();
    }


    uint32_t multiple_core_pair_map_index::get_k(uint32_t h) {

        return left_map->count(h) ? left_map->at(h) : 0;
    }

    uint32_t multiple_core_pair_map_index::get_h(uint32_t k) {
        return right_map->count(k) ? right_map->at(k) : 0;
    }

    void multiple_core_pair_map_index::merge_insert(const shared_ptr<scnu::multiple_core_pair_map_index> &other_index) {
        for(const auto &[k,h]:*other_index->get_right_map()){
            if(!right_map->count(k)){
                right_map->insert({k, h});
            }else if(right_map->at(k) < h){
                right_map->at(k) = h;
            }
        }
        for(const auto &[h, k]:*other_index->get_left_map()){
            if(!left_map->count(h)){
                left_map->insert({h, k});
            }else if(left_map->at(h) < k){
                left_map->at(h) = k;
            }
        }
    }

    void multiple_core_pair_map_index::merge_remove(const shared_ptr<scnu::multiple_core_pair_map_index> &other_index) {
        for(const auto &[k, h]:*other_index->get_right_map()){
            if(k == h){
                right_map->erase(k);
            }else{
                right_map->at(k) = h - 1;
                if(right_map->at(k) == 0){
                    right_map->erase(k);
                }
            }
        }
        for(const auto &[h, k]:*other_index->get_left_map()){
            if(k == h){
                left_map->erase(h);
            }else{
                left_map->at(h) = k - 1;
                if(left_map->at(h) == 0){
                    left_map->erase(h);
                }
            }
        }
    }

    void multiple_core_pair_map_index::insert(uint32_t k, uint32_t h) {
        if(k > h){
            left_insert(h, k);
        }else if (k < h){
            right_insert(k, h);
        }else{
            left_insert(h, k);
            right_insert(k, h);
        }
    }

    void multiple_core_pair_map_index::left_insert(uint32_t key, uint32_t value) {
        if (!left_map->count(key)) {
            left_map->insert({key,value});
        } else if (left_map->at(key) < value) {
            left_map->at(key) = value;
        }
    }

    void multiple_core_pair_map_index::right_insert(uint32_t key, uint32_t value) {
        if (!right_map->count(key)) {
            right_map->insert({key, value});
        } else if (right_map->at(key) < value) {
            right_map->at(key) = value;
        }
    }


    void multiple_core_pair_map_index::remove(uint32_t k, uint32_t h) {
        if(k > h){
            left_remove(h, k);
        }
        else if(k < h){
            right_remove(k, h);
        }else{
            left_remove(h, k);
            right_remove(k, h);
        }
    }

    void multiple_core_pair_map_index::left_remove(uint32_t key) {
        if (left_map->count(key)) {
            left_map->erase(key);
        }
    }

    void multiple_core_pair_map_index::left_remove(uint32_t key, uint32_t value) {
        if (!left_map->count(key)) {
            left_map->insert({key, value});
        } else if (left_map->at(key) > value) {
            left_map->at(key) = value;
        }
    }

    void multiple_core_pair_map_index::right_remove(uint32_t key) {
        if (right_map->count(key)) {
            right_map->erase(key);
        }
    }

    void multiple_core_pair_map_index::right_remove(uint32_t key, uint32_t value) {
        if (!right_map->count(key)) {
            right_map->insert({key, value});
        } else if (right_map->at(key) > value) {
            right_map->at(key) = value;
        }
    }

    void multiple_core_pair_map_index::set(uint32_t k, uint32_t h) {
        if(k > h){
            set_left(h, k);
        }
        else if(k < h){
            set_right(k, h);
        }else{
            set_left(h, k);
            set_right(k, h);
        }
    }

    void multiple_core_pair_map_index::set_left(uint32_t key, uint32_t value) {
        if(left_map->count(key)){
            left_map->at(key) = value;
        }
    }

    void multiple_core_pair_map_index::set_right(uint32_t key, uint32_t value) {
        if(right_map->count(key)){
            right_map->at(key) = value;
        }
    }
}