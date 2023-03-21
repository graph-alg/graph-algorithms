
#include "multiple_core/hierarchy_multiple_core_maintenance.h"

namespace scnu {

    void
    hierarchy_multiple_core_maintenance::assign(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                uint32_t value,
                                                const shared_ptr<scnu::thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto &iter = sub_begin; iter!=sub_end;++iter){
                    auto &u = *iter;
                    vertex_degree_map->at(u) = value;
                }
            });
        }
        pool->barrier();
    }


    void hierarchy_multiple_core_maintenance::insertion_assign(const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                               uint32_t k, uint32_t h, const shared_ptr<scnu::thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(candidate_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto &iter = sub_begin; iter!=sub_end;++iter){
                    auto &[u, u_degree] = *iter;
                    vertex_index_map->at(u)->insert(k, h);
                }
            });
        }
        pool->barrier();
    }

    void hierarchy_multiple_core_maintenance::removal_assign(const shared_ptr<unordered_set<uint32_t>> &removed_set,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                             uint32_t k, uint32_t h, const shared_ptr<scnu::thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(removed_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto &iter = sub_begin; iter!=sub_end;++iter){
                    auto &u = *iter;
                    vertex_index_map->at(u)->remove(k, h);
                }
            });
        }
        pool->barrier();
    }

    void hierarchy_multiple_core_maintenance::insertion_merge(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                              const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(new_vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto &iter = sub_begin; iter!=sub_end;++iter){
                    auto &[u, u_index] = *iter;
                    vertex_index_map->at(u)->merge_insert(u_index);
                    u_index->clear();
                }
            });
        }
        pool->barrier();
    }

    void hierarchy_multiple_core_maintenance::removal_merge(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                            const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(new_vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto &iter = sub_begin; iter!=sub_end;++iter){
                    auto &[u, u_index] = *iter;
                    vertex_index_map->at(u)->merge_remove(u_index);
                    u_index->clear();
                }
            });
        }
        pool->barrier();
    }


    void hierarchy_multiple_core_maintenance::init(const shared_ptr<temporal_graph>& G,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                   const shared_ptr<thread_pool>& pool){
        for(const auto &[u, u_vertex]:*G->get_vertex_map()){
            vertex_mutex_map->insert({u, make_shared<mutex>()});
            vertex_degree_map->insert({u, 0});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_mutex_map);
        for(uint32_t i = 0; i< thread_number;++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter!=sub_end;++iter){
                    auto u = iter->first;
                    vertex_mutex_map->at(u) = make_shared<mutex>();
                }
            });
        }
        pool->barrier();
    }

    uint32_t hierarchy_multiple_core_maintenance::compute_left_core_degree(const shared_ptr<temporal_graph> &G,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                           uint32_t u,
                                                                           uint32_t k,
                                                                           uint32_t h) {
        uint32_t degree = 0;
        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
            if (vertex_index_map->at(v)->count(k, h) && G->get_vertex(v)->get_neighbor_size() >= k + 1 && e_set->size() >= h) {
                ++degree;
            }
        }
        return degree;
    }

    uint32_t hierarchy_multiple_core_maintenance::compute_right_core_degree(const shared_ptr<temporal_graph> &G,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                           uint32_t u,
                                                                           uint32_t k,
                                                                           uint32_t h) {
        uint32_t degree = 0;
        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
            if (vertex_index_map->at(v)->count(k, h) && e_set->size() >= h + 1) {
                ++degree;
            }
        }
        return degree;
    }

    uint32_t hierarchy_multiple_core_maintenance::compute_core_degree(const shared_ptr<temporal_graph> &G,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                      uint32_t u,
                                                                      uint32_t k,
                                                                      uint32_t h) {
        uint32_t degree = 0;
        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
            if (vertex_index_map->at(v)->count(k, h) && e_set->size() >= h) {
                ++degree;
            }
        }
        return degree;
    }

    void hierarchy_multiple_core_maintenance::find_left_quasi_core(const shared_ptr<temporal_graph> &G,
                                                                   const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                   uint32_t k,
                                                                   uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        for (const auto &[u, u_degree]: *vertex_degree_map) {
            if(u_degree < k){
                vertex_set->insert(u);
            }
        }

        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                vertex_degree_map->erase(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*vertex_set){
                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_degree_map->count(v) && v_edge_set->size() >= h && vertex_degree_map->at(v) >= k) {
                        --vertex_degree_map->at(v);
                        if (vertex_degree_map->at(v) < k) {
                            next_evicted_set->insert(v);
                        }
                    }
                }
            }
            swap(*vertex_set, *next_evicted_set);
        }

        for(auto iter = edge_set->begin(); iter!=edge_set->end();){
            auto &e = *iter;
            ++iter;
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            if(!vertex_degree_map->count(u) || !vertex_degree_map->count(v)){
                edge_set->erase(e);
            }
        }
    }

    void hierarchy_multiple_core_maintenance::find_left_quasi_core(const shared_ptr<temporal_graph> &G,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                   const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                   uint32_t k,
                                                                   uint32_t h,
                                                                   const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(vertex_degree_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &[u, u_degree] = *iter;
                        if(u_degree < k){
                            sub_vertex_set->insert(u);
                        }
                    }
                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                vertex_degree_map->erase(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_removed_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &u = *iter;
                        for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (vertex_degree_map->count(v)  && v_edge_set->size() >= h && vertex_degree_map->at(v) >= k) {
                                vertex_mutex_map->at(v)->lock();
                                --vertex_degree_map->at(v);
                                vertex_mutex_map->at(v)->unlock();

                                if (vertex_degree_map->at(v) < k) {
                                    sub_next_evicted_set->insert(v);
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    next_evicted_set->merge(*sub_next_evicted_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            swap(*vertex_set, *next_evicted_set);
        }

        for(auto iter = edge_set->begin(); iter!=edge_set->end();){
            auto &e = *iter;
            ++iter;
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            if(!vertex_degree_map->count(u) || !vertex_degree_map->count(v)){
                edge_set->erase(e);
            }
        }
    }

    void hierarchy_multiple_core_maintenance::find_right_quasi_core(const shared_ptr<temporal_graph> &G,
                                                                    const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                    uint32_t k,
                                                                    uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        for (const auto &[u, u_degree]: *vertex_degree_map) {
            for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                if (vertex_degree_map->count(v) && e_set->size() == h - 1) {
                    --vertex_degree_map->at(u);
                }
            }
            if (vertex_degree_map->at(u) < k) {
                vertex_set->insert(u);
            }
        }

        for (const auto &u: *vertex_set) {
            vertex_degree_map->erase(u);
        }

        for(auto iter = edge_set->begin(); iter!=edge_set->end();){
            auto &e = *iter;
            ++iter;
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            if(!vertex_degree_map->count(u) || !vertex_degree_map->count(v)){
                edge_set->erase(e);
            }
        }
    }


    void hierarchy_multiple_core_maintenance::find_right_quasi_core(const shared_ptr<temporal_graph> &G,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                   const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                   uint32_t k,
                                                                   uint32_t h,
                                                                   const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(vertex_degree_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &[u, u_degree] = *iter;
                        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (vertex_degree_map->count(v) && e_set->size() == h - 1) {
                                --vertex_degree_map->at(u);
                            }
                        }
                        if (vertex_degree_map->at(u) < k) {
                            sub_vertex_set->insert(u);
                        }
                    }
                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        for (const auto &u: *vertex_set) {
            vertex_degree_map->erase(u);
        }

        for (auto iter = edge_set->begin(); iter != edge_set->end();) {
            auto &e = *iter;
            ++iter;
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            if (!vertex_degree_map->count(u) || !vertex_degree_map->count(v)) {
                edge_set->erase(e);
            }
        }
    }

    void hierarchy_multiple_core_maintenance::find_quasi_cores(const shared_ptr<temporal_graph> &G,
                                                               const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                               const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_vertex_degree_map,
                                                               uint32_t k,
                                                               uint32_t h,
                                                               const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>> &quasi_edge_map) {
        auto vertex_degree_map = quasi_vertex_degree_map->at({k, h});
        quasi_vertex_degree_map->erase({k, h});

        if (k == 1) {
            auto sub_edge_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>(*edge_set);
            auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
            find_right_quasi_core(G, sub_edge_set, sub_vertex_degree_map, k, h + 1);
            if(!sub_edge_set->empty()){
                quasi_edge_map->insert({{k, h + 1}, sub_edge_set});
                quasi_vertex_degree_map->insert({{k, h + 1}, sub_vertex_degree_map});
            }
        }

        auto sub_edge_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>(*edge_set);
        auto &sub_vertex_degree_map = vertex_degree_map;
        find_left_quasi_core(G, sub_edge_set, sub_vertex_degree_map, k + 1, h);

        if(!sub_edge_set->empty()){
            quasi_edge_map->insert({{k + 1, h}, sub_edge_set});
            quasi_vertex_degree_map->insert({{k + 1, h}, sub_vertex_degree_map});
        }
    }

    void hierarchy_multiple_core_maintenance::find_quasi_cores(const shared_ptr<temporal_graph> &G,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                               const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                               const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_vertex_degree_map,
                                                               uint32_t k,
                                                               uint32_t h,
                                                               const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>> &quasi_edge_map,
                                                               const shared_ptr<thread_pool>& pool) {
        auto vertex_degree_map = quasi_vertex_degree_map->at({k, h});
        quasi_vertex_degree_map->erase({k, h});

        if (k == 1) {
            auto sub_edge_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>(*edge_set);
            auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
            find_right_quasi_core(G, vertex_mutex_map, sub_edge_set, sub_vertex_degree_map, k, h + 1, pool);
            if(!sub_edge_set->empty()){
                quasi_edge_map->insert({{k, h + 1}, sub_edge_set});
                quasi_vertex_degree_map->insert({{k, h + 1}, sub_vertex_degree_map});
            }
        }

        auto sub_edge_set = make_shared<unordered_set<shared_ptr<temporal_edge>>>(*edge_set);
        auto &sub_vertex_degree_map = vertex_degree_map;
        find_left_quasi_core(G, vertex_mutex_map, sub_edge_set, sub_vertex_degree_map, k + 1, h, pool);

        if(!sub_edge_set->empty()){
            quasi_edge_map->insert({{k + 1, h}, sub_edge_set});
            quasi_vertex_degree_map->insert({{k + 1, h}, sub_vertex_degree_map});
        }
    }

    bool hierarchy_multiple_core_maintenance::left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                   const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                   uint32_t k,
                                                                   uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t count = 0;
        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            if (vertex_index_map->at(u)->count(k - 1, h) && vertex_index_map->at(v)->count(k - 1, h)) {
                if(!vertex_index_map->at(u)->count(k, h)){
                    if (!vertex_degree_map->count(u)) {
                        vertex_degree_map->insert({u, compute_left_core_degree(G, vertex_index_map, u, k - 1, h)});
                    }

                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }
                }

                if (!vertex_degree_map->count(v)) {
                    vertex_degree_map->insert({v, compute_left_core_degree(G, vertex_index_map, v, k - 1, h)});
                }

                if (vertex_degree_map->at(v) >= k) {
                    vertex_set->insert(v);
                }
            } else {
                ++count;
            }
        }

        if (count == edge_set->size()) {
            return false;
        }

        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto u = *vertex_set->begin();
            vertex_set->erase(u);

            uint32_t u_degree = 0;
            auto u_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                if (vertex_index_map->at(v)->count(k - 1, h) && G->get_vertex(v)->get_neighbor_size()>= k && e_set->size() >= h && !evicted_set->count(v)) {
                    if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v)) {
                        ++u_degree;
                    } else
                    {
                        if (!vertex_degree_map->count(v)) {
                            vertex_degree_map->insert({v, compute_left_core_degree(G, vertex_index_map, v, k - 1, h)});
                        }
                        if (vertex_degree_map->at(v) >= k) {
                            ++u_degree;
                            u_set->insert(v);
                        }
                    }
                }
            }
            if (u_degree >= k) {
                candidate_map->insert({u, u_degree});
                vertex_set->merge(*u_set);
            } else {
                remove_unsatisfied_vertices(G, u, candidate_map, evicted_set, k, h);
            }
        }
        return true;
    }

    bool hierarchy_multiple_core_maintenance::left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                   const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                   uint32_t k,
                                                                   uint32_t h,
                                                                   const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto count = make_shared<uint32_t>(0);
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        {

            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0;i < thread_number; ++i){
                pool->submit_task([=]{
                   auto &sub_begin = *location_vector->at(i);
                   auto &sub_end = *location_vector->at(i + 1);
                   auto sub_count = 0;
                   auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                   auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                   for(auto iter = sub_begin; iter!=sub_end; ++iter) {
                       auto &e = *iter;
                       auto u = e->get_source_vertex_id();
                       auto v = e->get_destination_vertex_id();

                       if (vertex_index_map->at(u)->count(k - 1, h) && vertex_index_map->at(v)->count(k - 1, h)) {
                           if (!vertex_index_map->at(u)->count(k, h)) {
                               vertex_mutex_map->at(u)->lock();
                               if (vertex_degree_map->at(u) == 0) {
                                   sub_clear_set->insert(u);
                                   vertex_degree_map->at(u) = compute_left_core_degree(G, vertex_index_map, u, k - 1,
                                                                                       h);
                               }
                               vertex_mutex_map->at(u)->unlock();

                               if (vertex_degree_map->at(u) >= k) {
                                   sub_vertex_set->insert(u);
                               }
                           }

                           vertex_mutex_map->at(v)->lock();
                           if (vertex_degree_map->at(v) == 0) {
                               sub_clear_set->insert(v);
                               vertex_degree_map->at(v) = compute_left_core_degree(G, vertex_index_map, v, k - 1, h);
                           }
                           vertex_mutex_map->at(v)->unlock();

                           if (vertex_degree_map->at(v) >= k) {
                               sub_vertex_set->insert(v);
                           }
                       } else {
                           ++sub_count;
                       }
                   }

                   global_mutex->lock();
                   *count += sub_count;
                   vertex_set->merge(*sub_vertex_set);
                   clear_set->merge(*sub_clear_set);
                   global_mutex->unlock();
                });
            }
            pool->barrier();
        }


        if (*count == edge_set->size()) {
            return false;
        }

        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto invalid_set = make_shared<unordered_set<uint32_t>>();
        auto visited_set = make_shared<unordered_set<uint32_t>>();

        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                   auto &sub_begin = *location_vector->at(i);
                   auto &sub_end = *location_vector->at(i + 1);

                   auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                   auto sub_valid_set = make_shared<unordered_set<uint32_t>>();
                   auto sub_layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                   for(auto iter = sub_begin; iter != sub_end; ++iter){
                       auto &u = *iter;
                       uint32_t u_degree = 0;
                       auto u_set = make_shared<unordered_set<uint32_t>>();
                       for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                           if (vertex_index_map->at(v)->count(k - 1, h) && G->get_vertex(v)->get_neighbor_size()>= k && e_set->size() >= h && !evicted_set->count(v)) {
                               if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v) || sub_layer_candidate_map->count(v)) {
                                   ++u_degree;
                               } else
                               {
                                   vertex_mutex_map->at(v)->lock();
                                   if (vertex_degree_map->at(v) == 0) {
                                       sub_clear_set->insert(v);
                                       vertex_degree_map->at(v) = compute_left_core_degree(G, vertex_index_map, v, k - 1, h);
                                   }
                                   vertex_mutex_map->at(v)->unlock();

                                   if (vertex_degree_map->at(v) >= k) {
                                       ++u_degree;
                                       if(!visited_set->count(v) && !vertex_set->count(v)){
                                           u_set->insert(v);
                                       }
                                   }
                               }
                           }
                       }
                       if (u_degree >= k) {
                           sub_layer_candidate_map->insert({u, u_degree});
                           sub_next_vertex_set->merge(*u_set);
                       } else {
                           sub_valid_set->insert(u);
                       }
                   }

                   global_mutex->lock();
                   layer_candidate_map->merge(*sub_layer_candidate_map);
                   next_vertex_set->merge(*sub_next_vertex_set);
                   invalid_set->merge(*sub_valid_set);
                   clear_set->merge(*sub_clear_set);
                   global_mutex->unlock();
                });
            }
            pool->barrier();
            candidate_map->merge(*layer_candidate_map);
            remove_unsatisfied_vertices(G, vertex_mutex_map, invalid_set, candidate_map, evicted_set,  k, h, pool);
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }

    bool hierarchy_multiple_core_maintenance::right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                    const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                    uint32_t k,
                                                                    uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        uint32_t count = 0;
        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            if (vertex_index_map->at(u)->count(k, h-1) && vertex_index_map->at(v)->count(k, h-1)) {
                if(G->get_edge_set(u, v) && G->get_edge_set(u, v)->size() >= h){
                    if (!vertex_degree_map->count(u)) {
                        vertex_degree_map->insert({u, compute_right_core_degree(G, vertex_index_map, u, k, h - 1)});
                    }

                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }

                    if (!vertex_degree_map->count(v)) {
                        vertex_degree_map->insert({v, compute_right_core_degree(G, vertex_index_map, v, k, h - 1)});
                    }

                    if (vertex_degree_map->at(v) >= k) {
                        vertex_set->insert(v);
                    }
                }
            } else {
                ++count;
            }
        }

        if (vertex_set->empty() || count == edge_set->size()) {
            return false;
        }

        for(const auto &u:*vertex_set){
            candidate_map->insert({u, UINT32_MAX});
        }
        return true;
    }

    bool hierarchy_multiple_core_maintenance::right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                    const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                    uint32_t k,
                                                                    uint32_t h,
                                                                    const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto count = make_shared<uint32_t>(0);
        auto location_vector = pool->split_task(edge_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
               auto &sub_begin = *location_vector->at(i);
               auto &sub_end = *location_vector->at(i + 1);

               auto sub_count = 0;
               auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
               auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

               for(auto iter = sub_begin; iter!=sub_end; ++iter){
                   auto &e = *iter;

                   auto u = e->get_source_vertex_id();
                   auto v = e->get_destination_vertex_id();

                   if (vertex_index_map->at(u)->count(k, h-1) && vertex_index_map->at(v)->count(k, h-1)) {
                       if(G->get_edge_set(u, v) && G->get_edge_set(u, v)->size() >= h){
                           vertex_mutex_map->at(u)->lock();
                           if (vertex_degree_map->at(u) == 0) {
                               sub_clear_set->insert(u);
                               vertex_degree_map->at(u) = compute_right_core_degree(G, vertex_index_map, u, k, h - 1);
                           }
                           vertex_mutex_map->at(u)->unlock();

                           if (vertex_degree_map->at(u) >= k) {
                               sub_vertex_set->insert(u);
                           }

                           vertex_mutex_map->at(v)->lock();
                           if (vertex_degree_map->at(v) == 0) {
                               sub_clear_set->insert(v);
                               vertex_degree_map->at(v) = compute_right_core_degree(G, vertex_index_map, v, k, h - 1);
                           }
                           vertex_mutex_map->at(v)->unlock();

                           if (vertex_degree_map->at(v) >= k) {
                               sub_vertex_set->insert(v);
                           }
                       }
                   } else {
                       ++sub_count;
                   }
               }

               global_mutex->lock();
               *count += sub_count;
               vertex_set->merge(*sub_vertex_set);
               clear_set->merge(*sub_clear_set);
               global_mutex->unlock();
            });
        }
        pool->barrier();

        if (*count == edge_set->size()) {
            assign(clear_set, vertex_degree_map, 0, pool);
            return false;
        }

        for(const auto &u:*vertex_set){
            candidate_map->insert({u, UINT32_MAX});
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }

    void hierarchy_multiple_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                                          uint32_t w,
                                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                          const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                                          uint32_t k,
                                                                          uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(w);

        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                candidate_map->erase(u);
                evicted_set->insert(u);
            }
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*vertex_set){
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (candidate_map->count(v) && e_set->size() >= h && candidate_map->at(v) >= k) {
                        --candidate_map->at(v);
                        if (candidate_map->at(v) < k) {
                            next_vertex_set->insert(v);
                        }
                    }
                }
            }
            swap(*vertex_set, *next_vertex_set);
        }
    }

    void hierarchy_multiple_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                          const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                          const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                                          uint32_t k,
                                                                          uint32_t h,
                                                                          const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                candidate_map->erase(u);
                evicted_set->insert(u);
            }
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &u = *iter;

                        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (candidate_map->count(v) && e_set->size() >= h && candidate_map->at(v) >= k) {

                                vertex_mutex_map->at(v)->lock();
                                --candidate_map->at(v);
                                vertex_mutex_map->at(v)->unlock();

                                if (candidate_map->at(v) < k) {
                                    sub_next_vertex_set->insert(v);
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    next_vertex_set->merge(*sub_next_vertex_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            swap(*vertex_set, *next_vertex_set);
        }
    }




    void hierarchy_multiple_core_maintenance::insert(const shared_ptr<temporal_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map) {
        auto new_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
        auto quasi_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        {
            auto affected_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                affected_set->insert(u);
                affected_set->insert(v);
                G->insert_edge(e);
            }

            quasi_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            auto quasi_degree_map = quasi_vertex_degree_map->at({1, 1});
            for (const auto &u: *affected_set) {
                quasi_degree_map->insert({u, G->get_vertex(u)->get_neighbor_size()});

                if(!vertex_index_map->count(u)){
                    vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    vertex_index_map->at(u)->insert(1,1);
                }
            }
            find_quasi_cores(G, edge_set, quasi_vertex_degree_map, 1, 1, quasi_edge_map);
        }

        while (!quasi_edge_map->empty()) {
            auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            for (const auto &p: *quasi_edge_map) {
                auto [k, h] = p.first;
                if (k == 1) {
                    auto sub_edge_set = quasi_edge_map->at({k, h});
                    auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto candidate_map = make_shared<unordered_map<uint32_t,uint32_t>>();
                    auto flag = right_candidate_graph(G, sub_edge_set, vertex_index_map, vertex_degree_map, candidate_map, k, h);

                    if (flag) {
                        for (const auto &[u, u_degree]: *candidate_map) {
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                            new_vertex_index_map->at(u)->insert(k, h);
                        }
                    }else{
                        flag_set->insert({k, h});
                    }

                } else {
                    auto sub_edge_set = quasi_edge_map->at({k, h});
                    auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto flag = left_candidate_graph(G, sub_edge_set, vertex_index_map, vertex_degree_map, candidate_map, k, h);

                    if (flag) {
                        for (const auto &[u, u_degree]: *candidate_map) {
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                            new_vertex_index_map->at(u)->insert(k, h);
                        }
                    }else{
                        flag_set->insert({k, h});
                    }
                }
            }
            for (const auto &[u, u_index]: *new_vertex_index_map) {
                vertex_index_map->at(u)->merge_insert(u_index);
                u_index->clear();
            }
            for(const auto &p:*quasi_edge_map){
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;
                if(!flag_set->count(p.first)){
                    find_quasi_cores(G, sub_edge_set, quasi_vertex_degree_map, k, h, next_quasi_edge_map);
                }
            }
            quasi_edge_map->clear();
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }
    }


    void hierarchy_multiple_core_maintenance::insert(const shared_ptr<temporal_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                     const shared_ptr<thread_pool>& pool) {
        auto new_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
        auto quasi_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t,uint32_t>>, hash_pair, equal_pair>>();
        {
            auto affected_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                affected_set->insert(u);
                affected_set->insert(v);
                G->insert_edge(e);
            }

            quasi_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            auto quasi_degree_map = quasi_vertex_degree_map->at({1, 1});

            for (const auto &u: *affected_set) {
                quasi_degree_map->insert({u, G->get_vertex(u)->get_neighbor_size()});

                if(!vertex_index_map->count(u)){
                    vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    vertex_index_map->at(u)->insert(1,1);
                }
            }
            find_quasi_cores(G, edge_set, quasi_vertex_degree_map, 1, 1, quasi_edge_map);
        }

        auto global_mutex = make_shared<mutex>();
        while (!quasi_edge_map->empty()) {
            auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            if(quasi_edge_map->size() > 1){
                for (const auto &p: *quasi_edge_map) {
                    if (p.first.first == 1) {
                        pool->submit_task([=] {
                            auto [k, h] = p.first;
                            auto sub_edge_set = quasi_edge_map->at({k, h});
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto flag = right_candidate_graph(G, sub_edge_set, vertex_index_map, vertex_degree_map,
                                                              candidate_map, k, h);

                            global_mutex->lock();
                            if (flag) {
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(k, h);
                                }
                            } else {
                                flag_set->insert({k, h});
                            }
                            global_mutex->unlock();
                        });
                    } else {
                        pool->submit_task([=] {
                            auto [k, h] = p.first;
                            auto sub_edge_set = quasi_edge_map->at({k, h});
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto flag = left_candidate_graph(G, sub_edge_set, vertex_index_map, vertex_degree_map,
                                                             candidate_map, k, h);
                            global_mutex->lock();
                            if (flag) {
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(k, h);
                                }
                            } else {
                                flag_set->insert({k, h});
                            }
                            global_mutex->unlock();
                        });
                    }
                }
                pool->barrier();
            }else{
                for (const auto &p: *quasi_edge_map) {
                    if (p.first.first == 1) {
                            auto [k, h] = p.first;
                            auto sub_edge_set = quasi_edge_map->at({k, h});
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto flag = right_candidate_graph(G, sub_edge_set, vertex_index_map, vertex_degree_map,
                                                              candidate_map, k, h);

                            if (flag) {
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(k, h);
                                }
                            } else {
                                flag_set->insert({k, h});
                            }
                    } else {
                        auto [k, h] = p.first;
                        auto sub_edge_set = quasi_edge_map->at({k, h});
                        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto flag = left_candidate_graph(G, sub_edge_set, vertex_index_map, vertex_degree_map,
                                                         candidate_map, k, h);

                        if (flag) {
                            for (const auto &[u, u_degree]: *candidate_map) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->insert(k, h);
                            }
                        } else {
                            flag_set->insert({k, h});
                        }
                    }
                }
            }
            for(const auto &p:*quasi_edge_map){
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;
                if(!flag_set->count(p.first)){
                    find_quasi_cores(G, sub_edge_set, quasi_vertex_degree_map, k, h, next_quasi_edge_map);
                }
            }
            for (const auto &[u, u_index]: *new_vertex_index_map) {
                vertex_index_map->at(u)->merge_insert(u_index);
                u_index->clear();
            }
            quasi_edge_map->clear();
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }
    }

    void hierarchy_multiple_core_maintenance::insert(const shared_ptr<temporal_graph> &G,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                     const shared_ptr<thread_pool>& pool) {
        auto new_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
        auto quasi_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t,uint32_t>>, hash_pair, equal_pair>>();
        {
            auto affected_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                affected_set->insert(u);
                affected_set->insert(v);
                G->insert_edge(e);
            }

            quasi_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            auto quasi_degree_map = quasi_vertex_degree_map->at({1, 1});

            for (const auto &u: *affected_set) {
                quasi_degree_map->insert({u, G->get_vertex(u)->get_neighbor_size()});

                if(!vertex_index_map->count(u)){
                    vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    vertex_index_map->at(u)->insert(1,1);

                    vertex_mutex_map->insert({u, make_shared<mutex>()});
                    vertex_degree_map->insert({u, 0});
                }
            }
            find_quasi_cores(G, vertex_mutex_map, edge_set, quasi_vertex_degree_map, 1, 1, quasi_edge_map, pool);
        }

        while (!quasi_edge_map->empty()) {
            auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            for (const auto &p: *quasi_edge_map) {
                    auto [k, h] = p.first;
                    if (k == 1) {
                        auto sub_edge_set = quasi_edge_map->at({k, h});
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto flag = right_candidate_graph(G, vertex_mutex_map, sub_edge_set, vertex_index_map, vertex_degree_map, candidate_map, k, h, pool);

                        if (flag) {
                            for(const auto &[u, degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                            }
                            insertion_assign(candidate_map, new_vertex_index_map, k, h, pool);
                        }else{
                            flag_set->insert({k, h});
                        }
                    } else {
                        auto sub_edge_set = quasi_edge_map->at({k, h});
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto flag = left_candidate_graph(G, vertex_mutex_map, sub_edge_set, vertex_index_map, vertex_degree_map, candidate_map, k,  h, pool);

                        if (flag) {
                            for(const auto &[u, degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                            }
                            insertion_assign(candidate_map, new_vertex_index_map, k, h, pool);
                        }else{
                            flag_set->insert({k, h});
                        }
                    }
            }
            for(const auto &p:*quasi_edge_map){
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;
                if(!flag_set->count(p.first)){
                    find_quasi_cores(G, vertex_mutex_map, sub_edge_set, quasi_vertex_degree_map, k, h, next_quasi_edge_map,pool);
                }
            }

            insertion_merge(new_vertex_index_map, vertex_index_map, pool);

            quasi_edge_map->clear();
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }
    }


    void hierarchy_multiple_core_maintenance::remove(const shared_ptr<temporal_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map) {
        auto new_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
        auto quasi_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        {
            auto affected_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                affected_set->insert(u);
                affected_set->insert(v);
            }
            quasi_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            auto quasi_degree_map = quasi_vertex_degree_map->at({1, 1});
            for (const auto &u: *affected_set) {
                quasi_degree_map->insert({u, G->get_vertex(u)->get_neighbor_size()});
            }
            find_quasi_cores(G, edge_set, quasi_vertex_degree_map, 1, 1, quasi_edge_map);
        }

        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        G->remove_edge_collection(edge_set, isolated_vertex_set);

        while (!quasi_edge_map->empty()) {
            auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            for (const auto &p: *quasi_edge_map) {
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;

                auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto removed_set = make_shared<unordered_set<uint32_t>>();
                auto flag = update_single_core(G, edge_set, vertex_index_map, vertex_degree_map, removed_set, k, h);

                if (flag) {
                    for (const auto &u: *removed_set) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                        }
                        new_vertex_index_map->at(u)->remove(k, h);
                    }
                }else{
                    flag_set->insert({k, h});
                }
            }
            G->insert_edge_collection(edge_set);
            for (const auto &p: *quasi_edge_map) {
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;
                if(!flag_set->count({k, h})){
                    find_quasi_cores(G, sub_edge_set, quasi_vertex_degree_map, k, h, next_quasi_edge_map);
                }
            }
            G->remove_edge_collection(edge_set);
            quasi_edge_map->clear();
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }

        for (const auto &[u, u_index]: *new_vertex_index_map) {
            vertex_index_map->at(u)->merge_remove(u_index);
            u_index->clear();
        }

        for (const auto &u: *isolated_vertex_set) {
            vertex_index_map->erase(u);
            new_vertex_index_map->erase(u);
        }
    }


    void hierarchy_multiple_core_maintenance::remove(const shared_ptr<temporal_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                     const shared_ptr<thread_pool>& pool) {
        auto new_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
        auto quasi_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        {
            auto affected_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                affected_set->insert(u);
                affected_set->insert(v);
            }
            quasi_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            auto quasi_degree_map = quasi_vertex_degree_map->at({1, 1});
            for (const auto &u: *affected_set) {
                quasi_degree_map->insert({u, G->get_vertex(u)->get_neighbor_size()});
            }
            find_quasi_cores(G, edge_set, quasi_vertex_degree_map, 1, 1, quasi_edge_map);
        }

        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        G->remove_edge_collection(edge_set, isolated_vertex_set);

        auto global_mutex = make_shared<mutex>();
        while (!quasi_edge_map->empty()) {
            auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            if(quasi_edge_map->size() > 1){
                for (const auto &p: *quasi_edge_map) {
                    pool->submit_task([=]{
                        auto [k, h] = p.first;
                        auto sub_edge_set = p.second;

                        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto removed_set = make_shared<unordered_set<uint32_t>>();
                        auto flag = update_single_core(G, edge_set, vertex_index_map,vertex_degree_map, removed_set, k, h);

                        global_mutex->lock();
                        if (flag) {
                            for (const auto &u: *removed_set) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(k, h);
                            }
                        }else{
                            flag_set->insert({k, h});
                        }
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
            }else{
                for (const auto &p: *quasi_edge_map) {
                    auto [k, h] = p.first;
                    auto sub_edge_set = p.second;

                    auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto removed_set = make_shared<unordered_set<uint32_t>>();
                    auto flag = update_single_core(G, edge_set, vertex_index_map, vertex_degree_map, removed_set, k, h);

                    if (flag) {
                        for (const auto &u: *removed_set) {
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                            new_vertex_index_map->at(u)->remove(k, h);
                        }
                    } else {
                        flag_set->insert({k, h});
                    }
                }
            }

            G->insert_edge_collection(edge_set);
            for (const auto &p: *quasi_edge_map) {
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;
                if(!flag_set->count({k, h})){
                    find_quasi_cores(G, sub_edge_set, quasi_vertex_degree_map, k, h, next_quasi_edge_map);
                }
            }
            G->remove_edge_collection(edge_set);
            quasi_edge_map->clear();
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }

        for (const auto &[u, u_index]: *new_vertex_index_map) {
            vertex_index_map->at(u)->merge_remove(u_index);
            u_index->clear();
        }

        for (const auto &u: *isolated_vertex_set) {
            vertex_index_map->erase(u);
        }
    }

    void hierarchy_multiple_core_maintenance::remove(const shared_ptr<temporal_graph> &G,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                     const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                     const shared_ptr<thread_pool>& pool) {
        auto new_vertex_index_map = make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
        auto quasi_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        {
            auto affected_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();
                affected_set->insert(u);
                affected_set->insert(v);
            }

            quasi_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            auto quasi_degree_map = quasi_vertex_degree_map->at({1, 1});

            for (const auto &u: *affected_set) {
                quasi_degree_map->insert({u, G->get_vertex(u)->get_neighbor_size()});
            }
            find_quasi_cores(G, vertex_mutex_map, edge_set, quasi_vertex_degree_map, 1, 1, quasi_edge_map, pool);
        }

        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        G->remove_edge_collection(edge_set, isolated_vertex_set);

        while (!quasi_edge_map->empty()) {
            auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<temporal_edge>>>, hash_pair, equal_pair>>();
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            for (const auto &p: *quasi_edge_map) {
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;

                auto removed_set = make_shared<unordered_set<uint32_t>>();
                auto flag = update_single_core(G, vertex_mutex_map, edge_set, vertex_index_map, vertex_degree_map, removed_set, k, h, pool);

                if (flag) {
                    for(const auto &u:*removed_set){
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                        }
                    }
                    removal_assign(removed_set, new_vertex_index_map, k, h, pool);
                }else{
                    flag_set->insert({k, h});
                }
            }
            G->insert_edge_collection(edge_set);
            for (const auto &p: *quasi_edge_map) {
                auto [k, h] = p.first;
                auto sub_edge_set = p.second;
                if(!flag_set->count({k, h})){
                    find_quasi_cores(G, vertex_mutex_map, sub_edge_set, quasi_vertex_degree_map, k, h, next_quasi_edge_map, pool);
                }
            }
            G->remove_edge_collection(edge_set);
            quasi_edge_map->clear();
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }

        removal_merge(new_vertex_index_map, vertex_index_map, pool);

        for (const auto &u: *isolated_vertex_set) {
            vertex_index_map->erase(u);
            vertex_mutex_map->erase(u);
            vertex_degree_map->erase(u);
        }
    }

    bool hierarchy_multiple_core_maintenance::update_single_core(const shared_ptr<temporal_graph> &G,
                                                                 const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                 const shared_ptr<unordered_set<uint32_t>> & removed_set,
                                                                 uint32_t k,
                                                                 uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t count = 0;
        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            if (vertex_index_map->at(u)->count(k, h) && vertex_index_map->at(v)->count(k, h)) {
                auto sub_edge_set = G->get_edge_set(u,v);
                if(!sub_edge_set || sub_edge_set->size() < h){
                    if(G->get_vertex(u)){
                        if (!vertex_degree_map->count(u)) {
                            vertex_degree_map->insert({u, compute_core_degree(G, vertex_index_map, u, k, h)});
                        }
                        if (vertex_degree_map->at(u) < k) {
                            vertex_set->insert(u);
                        }
                    }

                    if(G->get_vertex(v)){
                        if (!vertex_degree_map->count(v)) {
                            vertex_degree_map->insert({v, compute_core_degree(G, vertex_index_map, v, k, h)});
                        }
                        if (vertex_degree_map->at(v) < k) {
                            vertex_set->insert(v);
                        }
                    }
                }
            } else {
                ++count;
            }
        }

        if (count == edge_set->size()) {
            return false;
        }

        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for(const auto&u:*vertex_set){
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_index_map->at(v)->count(k, h) && e_set->size() >= h && !vertex_set->count(v) &&!removed_set->count(v)) {
                        if (!vertex_degree_map->count(v)) {
                            vertex_degree_map->insert({v, compute_core_degree(G, vertex_index_map, v, k, h)});
                        }

                        --vertex_degree_map->at(v);
                        if (vertex_degree_map->at(v) < k) {
                            next_vertex_set->insert(v);
                        }
                    }
                }
            }
            removed_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return true;
    }

    bool hierarchy_multiple_core_maintenance::update_single_core(const shared_ptr<temporal_graph> &G,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                 const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                 const shared_ptr<unordered_set<uint32_t>> & removed_set,
                                                                 uint32_t k,
                                                                 uint32_t h,
                                                                 const shared_ptr<thread_pool>& pool) {

        uint32_t const thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto count = make_shared<uint32_t>(0);
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_count = 0;
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto sub_begin = *location_vector->at(i);
                    auto sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &e = *iter;
                        auto u = e->get_source_vertex_id();
                        auto v = e->get_destination_vertex_id();

                        if (vertex_index_map->at(u)->count(k, h) && vertex_index_map->at(v)->count(k, h)) {
                            auto sub_edge_set = G->get_edge_set(u, v);
                            if (!sub_edge_set || sub_edge_set->size() < h) {
                                if (G->get_vertex(u)) {

                                    vertex_mutex_map->at(u)->lock();
                                    if (vertex_degree_map->at(u) == 0) {
                                        sub_clear_set->insert(u);
                                        vertex_degree_map->at(u) = compute_core_degree(G, vertex_index_map, u, k, h);
                                    }
                                    vertex_mutex_map->at(u)->unlock();

                                    if (vertex_degree_map->at(u) < k) {
                                        sub_vertex_set->insert(u);
                                    }
                                }

                                if (G->get_vertex(v)) {

                                    vertex_mutex_map->at(v)->lock();
                                    if (vertex_degree_map->at(v) == 0) {
                                        sub_clear_set->insert(v);
                                        vertex_degree_map->at(v) = compute_core_degree(G, vertex_index_map, v, k, h);
                                    }
                                    vertex_mutex_map->at(v)->unlock();

                                    if (vertex_degree_map->at(v) < k) {
                                        sub_vertex_set->insert(v);
                                    }
                                }
                            }
                        } else {
                            ++sub_count;
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    *count += sub_count;
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        if (*count == edge_set->size()) {
            return false;
        }

        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto sub_begin = *location_vector->at(i);
                    auto sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;

                        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (vertex_index_map->at(v)->count(k, h) && e_set->size() >= h && !vertex_set->count(v) && !removed_set->count(v)) {

                                vertex_mutex_map->at(v)->lock();
                                if (vertex_degree_map->at(v) == 0) {
                                    sub_clear_set->insert(v);
                                    vertex_degree_map->at(v) = compute_core_degree(G, vertex_index_map, v, k, h);
                                }
                                --vertex_degree_map->at(v);
                                vertex_mutex_map->at(v)->unlock();

                                if (vertex_degree_map->at(v) < k) {
                                    sub_next_vertex_set->insert(v);
                                }
                            }
                        }
                    }

                    global_mutex->lock();
                    next_vertex_set->merge(*sub_next_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            removed_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }
}