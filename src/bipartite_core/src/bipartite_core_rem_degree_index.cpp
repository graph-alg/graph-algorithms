
#include "bipartite_core/bipartite_core_rem_degree_index.h"

namespace scnu{
    bipartite_core_rem_degree_index::bipartite_core_rem_degree_index():
            middle_map(make_shared<unordered_map<uint32_t, uint32_t>>()),
            left_map(make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>()),
            right_map(make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>()) {

    }

    bipartite_core_rem_degree_index::bipartite_core_rem_degree_index(const shared_ptr<scnu::bipartite_core_rem_degree_index> &other_core_rem_degree_index)
            :bipartite_core_rem_degree_index(){
        middle_map = container_copy::to_unordered_map<uint32_t, uint32_t>(other_core_rem_degree_index->get_middle_map());
        left_map = container_copy::to_unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>(other_core_rem_degree_index->get_left_map());
        right_map = container_copy::to_unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>(other_core_rem_degree_index->get_right_map());

        for(const auto &[k, k_map]:*other_core_rem_degree_index->get_left_map()){
            left_map->at(k) = container_copy::to_unordered_map<uint32_t, uint32_t>(k_map);
        }

        for(const auto &[k, k_map]:*other_core_rem_degree_index->get_right_map()){
            right_map->at(k) = container_copy::to_unordered_map<uint32_t, uint32_t>(k_map);
        }
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>
    bipartite_core_rem_degree_index::get_left_map() {
        return left_map;
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>>
    bipartite_core_rem_degree_index::get_left_map(uint32_t k) {
        return left_map->at(k);
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>> bipartite_core_rem_degree_index::get_middle_map() {
        return middle_map;
    }

    shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>
    bipartite_core_rem_degree_index::get_right_map() {
        return right_map;
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>>
    bipartite_core_rem_degree_index::get_right_map(uint32_t k) {
        return right_map->at(k);
    }

    void bipartite_core_rem_degree_index::insert_left(uint32_t k,
                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &degree_map) {
        left_map->insert({k, degree_map});
    }

    void bipartite_core_rem_degree_index::insert_right(uint32_t k,
                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &degree_map) {
        right_map->insert({k, degree_map});
    }
}