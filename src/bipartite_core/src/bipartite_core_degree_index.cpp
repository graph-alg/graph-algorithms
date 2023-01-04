
#include "bipartite_core/bipartite_core_degree_index.h"

namespace scnu{
    bipartite_core_degree_index::bipartite_core_degree_index()
        :left_map(make_shared<unordered_map<uint32_t, uint32_t>>()),right_map(make_shared<unordered_map<uint32_t, uint32_t>>()){

    }

    bipartite_core_degree_index::bipartite_core_degree_index(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_left_map,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_right_map)
                                                             :left_map(other_left_map), right_map(other_right_map){

    }

    bipartite_core_degree_index::bipartite_core_degree_index(const shared_ptr<scnu::bipartite_core_degree_index> &other_index)
            :bipartite_core_degree_index(other_index->get_left_map(), other_index->get_right_map()){

    }

    bool bipartite_core_degree_index::count_left_vertex(uint32_t l) {
        return left_map->count(l);
    }

    bool bipartite_core_degree_index::count_right_vertex(uint32_t r) {
        return right_map->count(r);
    }

    bool bipartite_core_degree_index::empty() {
        return left_map->empty() || right_map->empty();
    }

    uint32_t bipartite_core_degree_index::get_left_degree(uint32_t l) {
        return left_map->count(l) ? left_map->at(l):0;
    }

    uint32_t bipartite_core_degree_index::get_right_degree(uint32_t r) {
        return right_map->count(r) ? right_map->at(r):0;
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>> bipartite_core_degree_index::get_left_map() {
        return left_map;
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>> bipartite_core_degree_index::get_right_map() {
        return right_map;
    }

    void bipartite_core_degree_index::insert_left_degree(uint32_t l, uint32_t degree) {
        left_map->insert({l, degree});
    }

    void bipartite_core_degree_index::insert_right_degree(uint32_t r, uint32_t degree) {
        right_map->insert({r, degree});
    }

    void bipartite_core_degree_index::insert_left_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_left_map) {
        copy(other_left_map->begin(), other_left_map->end(), inserter(*left_map, left_map->end()));
    }

    void bipartite_core_degree_index::insert_right_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_right_map) {
        copy(other_right_map->begin(), other_right_map->end(), inserter(*right_map, right_map->end()));
    }

    void bipartite_core_degree_index::update_left_degree(uint32_t l, uint32_t degree) {
        left_map->at(l) = degree;
    }

    void bipartite_core_degree_index::update_right_degree(uint32_t r, uint32_t degree) {
        right_map->at(r) = degree;
    }

    void bipartite_core_degree_index::update_left_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>> &other_left_map) {
        for(const auto &[l, degree]:*other_left_map){
            left_map->at(l) = degree;
        }
    }

    void bipartite_core_degree_index::update_right_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>> &other_right_map) {
        for(const auto &[r, degree]:*other_right_map){
            right_map->at(r) = degree;
        }
    }
}