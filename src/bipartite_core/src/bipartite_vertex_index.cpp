
#include "bipartite_core/bipartite_core_branch_store_index.h"

namespace scnu {
    bipartite_core_branch_store_index::bipartite_core_branch_store_index()
            : left_map(make_shared<unordered_map<uint32_t, uint32_t>>()),
              right_map(make_shared<unordered_map<uint32_t, uint32_t>>()) {
    }

    bipartite_core_branch_store_index::bipartite_core_branch_store_index(
            const shared_ptr<bipartite_core_branch_store_index> &other_bipartite_vertex_index)
            : bipartite_core_branch_store_index() {
        for (const auto &[j, i]: *other_bipartite_vertex_index->left_map) {
            left_map->insert({j, i});
        }
        for (const auto &[i, j]: *other_bipartite_vertex_index->right_map) {
            right_map->insert({i, j});
        }
    }

    void bipartite_core_branch_store_index::clear() {
        left_map->clear();
        right_map->clear();
    }

    bool bipartite_core_branch_store_index::compare(const shared_ptr<bipartite_core_branch_store_index> &other_multiple_core_map_index) {
        auto other_left_map = other_multiple_core_map_index->get_left_map();
        auto other_right_map = other_multiple_core_map_index->get_right_map();

        return container_compare::same_associative_map(left_map, other_left_map)
               && container_compare::same_associative_map(right_map, other_right_map);
    }

    bool bipartite_core_branch_store_index::count(uint32_t i, uint32_t j) {
        return i > j ? (left_map->count(j) && left_map->at(j) >= i) : (right_map->count(i) && right_map->at(i) >= j);
    }

    bool bipartite_core_branch_store_index::empty() {
        return left_map->empty() && right_map->empty();
    }

    bool bipartite_core_branch_store_index::equal(uint32_t i, uint32_t j) {
        if (i > j) {
            return left_map->count(j) && left_map->at(j) == i;
        } else if (i < j) {
            return right_map->count(i) && right_map->at(i) == j;
        } else {
            return (left_map->count(j) && left_map->at(j) == i) || (right_map->count(i) && right_map->at(i) == j);
        }
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>> bipartite_core_branch_store_index::get_left_map() {
        return left_map;
    }

    shared_ptr<unordered_map<uint32_t, uint32_t>> bipartite_core_branch_store_index::get_right_map() {
        return right_map;
    }

    uint32_t bipartite_core_branch_store_index::get_memory_cost() {
        return sizeof(*left_map->begin()) * left_map->size()
               + sizeof(*right_map->begin()) * right_map->size();
    }

    uint32_t bipartite_core_branch_store_index::get_k() {
        return left_map->size();
    }


    uint32_t bipartite_core_branch_store_index::get_i(uint32_t j) {

        return left_map->count(j) ? left_map->at(j) : 0;
    }

    uint32_t bipartite_core_branch_store_index::get_j(uint32_t i) {
        return right_map->count(i) ? right_map->at(i) : 0;
    }

    void bipartite_core_branch_store_index::merge_insert(const shared_ptr<scnu::bipartite_core_branch_store_index> &other_index) {
        for (const auto &[i, j]: *other_index->get_right_map()) {
            if (!right_map->count(i)) {
                right_map->insert({i, j});
            } else if (right_map->at(i) < j) {
                right_map->at(i) = j;
            }
        }
        for (const auto &[j, i]: *other_index->get_left_map()) {
            if (!left_map->count(j)) {
                left_map->insert({j, i});
            } else if (left_map->at(j) < i) {
                left_map->at(j) = i;
            }
        }
    }

    void bipartite_core_branch_store_index::merge_remove(const shared_ptr<scnu::bipartite_core_branch_store_index> &other_index) {
        for (const auto &[i, j]: *other_index->get_right_map()) {
            if (i == j) {
                right_map->erase(i);
            } else {
                right_map->at(i) = j - 1;
                if (right_map->at(i) == 0) {
                    right_map->erase(i);
                }
            }
        }
        for (const auto &[j, i]: *other_index->get_left_map()) {
            if (i == j) {
                left_map->erase(j);
            } else {
                left_map->at(j) = i - 1;
                if (left_map->at(j) == 0) {
                    left_map->erase(j);
                }
            }
        }
    }

    void bipartite_core_branch_store_index::insert(uint32_t i, uint32_t j) {
        if (i > j) {
            left_insert(j, i);
        } else if (i < j) {
            right_insert(i, j);
        } else {
            left_insert(j, i);
            right_insert(i, j);
        }
    }

    void bipartite_core_branch_store_index::left_insert(uint32_t key, uint32_t value) {
        if (!left_map->count(key)) {
            left_map->insert({key, value});
        } else if (left_map->at(key) < value) {
            left_map->at(key) = value;
        }
    }

    void bipartite_core_branch_store_index::right_insert(uint32_t key, uint32_t value) {
        if (!right_map->count(key)) {
            right_map->insert({key, value});
        } else if (right_map->at(key) < value) {
            right_map->at(key) = value;
        }
    }


    void bipartite_core_branch_store_index::remove(uint32_t i, uint32_t j) {
        if (i > j) {
            left_remove(j, i);
        } else if (i < j) {
            right_remove(i, j);
        } else {
            left_remove(j, i);
            right_remove(i, j);
        }
    }

    void bipartite_core_branch_store_index::left_remove(uint32_t key) {
        if (left_map->count(key)) {
            left_map->erase(key);
        }
    }

    void bipartite_core_branch_store_index::left_remove(uint32_t key, uint32_t value) {
        if (!left_map->count(key)) {
            left_map->insert({key, value});
        } else if (left_map->at(key) > value) {
            left_map->at(key) = value;
        }
    }

    void bipartite_core_branch_store_index::right_remove(uint32_t key) {
        if (right_map->count(key)) {
            right_map->erase(key);
        }
    }

    void bipartite_core_branch_store_index::right_remove(uint32_t key, uint32_t value) {
        if (!right_map->count(key)) {
            right_map->insert({key, value});
        } else if (right_map->at(key) > value) {
            right_map->at(key) = value;
        }
    }

    void bipartite_core_branch_store_index::set(uint32_t i, uint32_t j) {
        if (i > j) {
            set_left(j, i);
        } else if (i < j) {
            set_right(i, j);
        } else {
            set_left(j, i);
            set_right(i, j);
        }
    }

    void bipartite_core_branch_store_index::set_left(uint32_t key, uint32_t value) {
        if (left_map->count(key)) {
            left_map->at(key) = value;
        }
    }

    void bipartite_core_branch_store_index::set_right(uint32_t key, uint32_t value) {
        if (right_map->count(key)) {
            right_map->at(key) = value;
        }
    }
}