
#include "bipartite_core/bipartite_core_order_index.h"

namespace scnu{
    bipartite_core_order_index::bipartite_core_order_index() {
        middle_map = make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>();
        left_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>>>();
        right_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>>>();
    }

    bipartite_core_order_index::bipartite_core_order_index(const shared_ptr<bipartite_core_order_index>& other_core_order_index):
    bipartite_core_order_index(){
        for(const auto &[k, k_list]:* other_core_order_index->get_middle_map()){
            middle_map->insert({k, make_shared<extend_list<int, uint32_t>>(k_list)});
        }

        for(const auto &[k, k_map]:*other_core_order_index->get_left_map()){
            left_map->insert({k, make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>()});
            for(const auto &[i, i_list]:*k_map){
                left_map->at(k)->insert({i, make_shared<extend_list<int, uint32_t>>(i_list)});
            }
        }

        for(const auto &[k, k_map]:*other_core_order_index->get_right_map()){
            right_map->insert({k, make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>()});
            for(const auto &[j, j_list]:*k_map){
                right_map->at(k)->insert({j, make_shared<extend_list<int, uint32_t>>(j_list)});
            }
        }
    }


    shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>>>
    bipartite_core_order_index::get_left_map() {
        return left_map;
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>
    bipartite_core_order_index::get_left_map(uint32_t k) {
        return left_map->at(k);
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>
    bipartite_core_order_index::get_middle_map() {
        return middle_map;
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>>>
    bipartite_core_order_index::get_right_map() {
        return right_map;
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>
    bipartite_core_order_index::get_right_map(uint32_t k) {
        return right_map->at(k);
    }

    void bipartite_core_order_index::insert_left(uint32_t k,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &order_map) {
        left_map->insert({k, order_map});
    }

    void bipartite_core_order_index::insert_right(uint32_t k,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &order_map) {
        right_map->insert({k, order_map});
    }
}