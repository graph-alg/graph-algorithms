
#include "multiple_core/branch_multiple_core_maintenance.h"

namespace scnu {

    void branch_multiple_core_maintenance::assign(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
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


    void branch_multiple_core_maintenance::insertion_assign(const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                            uint32_t k,
                                                            uint32_t h,
                                                            const shared_ptr<scnu::thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();

        auto location_vector = pool->split_task(candidate_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto &iter = sub_begin; iter!=sub_end;++iter){
                    auto &[u, u_degree] = *iter;
                    if(!vertex_index_map->at(u)){
                        vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                    }
                    vertex_index_map->at(u)->insert(k, h);
                }
            });
        }

        pool->barrier();
    }

    void branch_multiple_core_maintenance::removal_assign(const shared_ptr<unordered_set<uint32_t>> &removed_set,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                          uint32_t k,
                                                          uint32_t h,
                                                          const shared_ptr<scnu::thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();

        auto location_vector = pool->split_task(removed_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto &iter = sub_begin; iter!=sub_end;++iter){
                    auto &u = *iter;
                    if(!vertex_index_map->at(u)){
                        vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                    }

                    vertex_index_map->at(u)->remove(k, h);
                }
            });
        }
        pool->barrier();
    }

    void branch_multiple_core_maintenance::insertion_merge(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
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

    void branch_multiple_core_maintenance::removal_merge(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
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

    uint32_t branch_multiple_core_maintenance::compute_left_core_degree(const shared_ptr<temporal_graph> &G,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                        uint32_t u,
                                                                        uint32_t k,
                                                                        uint32_t h) {
        uint32_t degree = 0;
        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
            if ((vertex_index_map->at(v)->count(k, h) || previous_candidate_map->count(v)) && e_set->size() >= h) {
                ++degree;
            }
        }
        return degree;
    }


    uint32_t branch_multiple_core_maintenance::compute_left_core_degree(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            uint32_t u,
            uint32_t k,
            uint32_t h) {
        uint32_t degree = 0;
        auto edge_vertex_map = vertex_edge_index_map->at(u);
        for (auto iter = edge_vertex_map->lower_bound(h); iter != edge_vertex_map->end(); ++iter) {
            for (const auto &v: *iter->second) {
                if (vertex_index_map->at(v)->count(k, h) || previous_candidate_map->count(v)) {
                    ++degree;
                }
            }
        }
        return degree;
    }

    uint32_t branch_multiple_core_maintenance::compute_middle_core_degree(const shared_ptr<temporal_graph> &G,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                          uint32_t u,
                                                                          uint32_t k,
                                                                          uint32_t h) {
        uint32_t degree = 0;
        auto neighbor_map = G->get_vertex(u)->get_neighbor_map();
        if (neighbor_map->size() < k + 1) {
            return degree;
        }
        for (const auto &[v, e_set]: *neighbor_map) {
            if ((vertex_index_map->at(v)->count(k, h) || previous_candidate_map->count(v)) && e_set->size() >= h + 1) {
                ++degree;
            }
        }
        return degree;
    }

    uint32_t branch_multiple_core_maintenance::compute_middle_core_degree(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            uint32_t u, uint32_t k, uint32_t h) {
        uint32_t degree = 0;
        auto edge_vertex_map = vertex_edge_index_map->at(u);
        for (auto iter = edge_vertex_map->lower_bound(h + 1); iter != edge_vertex_map->end(); ++iter) {
            for (const auto &v: *iter->second) {
                if (vertex_index_map->at(v)->count(k, h) || previous_candidate_map->count(v)) {
                    ++degree;
                }
            }
        }
        return degree;
    }

    uint32_t branch_multiple_core_maintenance::compute_right_core_degree(const shared_ptr<temporal_graph> &G,
                                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                         uint32_t u,
                                                                         uint32_t k,
                                                                         uint32_t h) {
        uint32_t degree = 0;
        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
            if ((vertex_index_map->at(v)->count(k, h) || previous_candidate_map->count(v)) && e_set->size() >= h + 1) {
                ++degree;
            }
        }
        return degree;
    }

    uint32_t branch_multiple_core_maintenance::compute_right_core_degree(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            uint32_t u,
            uint32_t k,
            uint32_t h) {
        uint32_t degree = 0;
        auto edge_vertex_map = vertex_edge_index_map->at(u);
        for (auto iter = edge_vertex_map->lower_bound(h + 1); iter != edge_vertex_map->end(); ++iter) {
            for (const auto &v: *iter->second) {
                if (vertex_index_map->at(v)->count(k, h) || previous_candidate_map->count(v)) {
                    ++degree;
                }
            }
        }
        return degree;
    }


    uint32_t branch_multiple_core_maintenance::compute_core_degree(const shared_ptr<temporal_graph> &G,
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

    uint32_t branch_multiple_core_maintenance::compute_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                   uint32_t u,
                                                                   uint32_t k,
                                                                   uint32_t h) {
        uint32_t degree = 0;
        auto edge_map = vertex_edge_index_map->at(u);
        for (auto iter = edge_map->lower_bound(h); iter != edge_map->end(); ++iter) {
            for (const auto &v: *iter->second)
                if (vertex_index_map->at(v)->count(k, h)) {
                    ++degree;
                }
        }
        return degree;
    }

    uint32_t branch_multiple_core_maintenance::find_max_delta(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map){
        uint32_t max_delta = 0;
        for (const auto &[u, u_index]: *vertex_index_map) {
            auto size = u_index->get_left_map()->size();
            if (size > max_delta) {
                max_delta = size;
            }
        }
        return  max_delta;
    }


    uint32_t branch_multiple_core_maintenance::find_max_delta(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                              const shared_ptr<thread_pool>& pool){
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto max_delta = make_shared<uint32_t>(0);
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=] {
                uint32_t sub_max_delta = 0;

                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &[u, u_index] = *iter;
                    auto delta = u_index->get_left_map()->size();
                    if(delta > sub_max_delta)
                    {
                        sub_max_delta = delta;
                    }
                }

                global_mutex->lock();
                if(sub_max_delta > *max_delta){
                    *max_delta = sub_max_delta;
                }
                global_mutex->unlock();
            });
        }
        pool->barrier();
        return  *max_delta;
    }

    uint32_t branch_multiple_core_maintenance::find_max_delta(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map){
        uint32_t max_delta = 0;
        for (const auto &u: *vertex_set) {
            auto delta = vertex_index_map->at(u)->get_left_map()->size();
            if (delta > max_delta) {
                max_delta = delta;
            }
        }
        return  max_delta;
    }

    uint32_t branch_multiple_core_maintenance::find_max_delta(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                              const shared_ptr<thread_pool>& pool){
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto max_delta = make_shared<uint32_t>(0);
        auto location_vector = pool->split_task(vertex_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=] {
                uint32_t sub_max_delta = 0;

                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &u = *iter;
                    auto u_index = vertex_index_map->at(u);
                    auto delta = u_index->get_left_map()->size();
                    if(delta > sub_max_delta)
                    {
                        sub_max_delta = delta;
                    }
                }

                global_mutex->lock();
                if(sub_max_delta > *max_delta){
                    *max_delta = sub_max_delta;
                }
                global_mutex->unlock();
            });
        }
        pool->barrier();
        return  *max_delta;
    }

    uint32_t branch_multiple_core_maintenance::find_max_k(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                          uint32_t h){
        uint32_t max_k = 0;
        for (const auto &u: *vertex_set) {
            auto u_index = vertex_index_map->at(u);
            auto k = u_index->get_k(h);
            if (k > max_k) {
                max_k = k;
            }
        }
        return max_k;
    }

    uint32_t branch_multiple_core_maintenance::find_max_k(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                              uint32_t h,
                                                              const shared_ptr<thread_pool>& pool){
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto max_k = make_shared<uint32_t>(0);
        auto location_vector = pool->split_task(vertex_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                uint32_t  sub_max_k = 0;

                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = * location_vector->at(i + 1);

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &u = *iter;
                    auto &u_index = vertex_index_map->at(u);
                    auto k = u_index->get_i(h);
                    if (k > sub_max_k) {
                        sub_max_k = k;
                    }
                }

                global_mutex->lock();
                if(sub_max_k > *max_k){
                    *max_k = sub_max_k;
                }
                global_mutex->unlock();
            });
        }
        pool->barrier();
        return *max_k;
    }

    uint32_t branch_multiple_core_maintenance::find_max_h(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                          uint32_t k){
        uint32_t max_h = 0;
        for (const auto &u: *vertex_set) {
            auto u_index = vertex_index_map->at(u);
            auto h = u_index->get_h(k);
            if (h > max_h) {
                max_h = h;
            }
        }
        return max_h;
    }

    uint32_t branch_multiple_core_maintenance::find_max_h(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                          uint32_t k,
                                                          const shared_ptr<thread_pool>& pool){
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto max_h = make_shared<uint32_t>(0);
        auto location_vector = pool->split_task(vertex_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                uint32_t  sub_max_h = 0;

                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = * location_vector->at(i + 1);

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &u = *iter;
                    auto &u_index = vertex_index_map->at(u);
                    auto h = u_index->get_j(k);
                    if (h > sub_max_h) {
                        sub_max_h = h;
                    }
                }

                global_mutex->lock();
                if(sub_max_h > *max_h){
                    *max_h = sub_max_h;
                }
                global_mutex->unlock();
            });
        }
        pool->barrier();
        return *max_h;
    }



    bool branch_multiple_core_maintenance::left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                const shared_ptr<unordered_map<uint32_t , uint32_t>> &previous_candidate_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                uint32_t k,
                                                                uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t count = 0;
        for (auto iter = affected_set->begin(); iter != affected_set->end();) {
            auto &u = *iter;
            ++iter;
            if (vertex_index_map->at(u)->count(k - 1, h) || previous_candidate_map->count(u)) {
                if (!vertex_index_map->at(u)->count(k, h)) {
                    vertex_degree_map->insert({u, compute_left_core_degree(G, vertex_index_map,
                                                                                 previous_candidate_map,
                                                                                 u, k - 1, h)});

                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }else {
                        ++count;
                    }
                }
            } else{
                ++count;
            }
        }
        if(count == affected_set->size()){
            return false;
        }


        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                uint32_t u_degree = 0;
                auto u_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if ((vertex_index_map->at(v)->count(k - 1, h) || previous_candidate_map->count(v)) &&
                        e_set->size() >= h && !evicted_set->count(v)) {
                        if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v)) {
                            ++u_degree;
                        } else {
                            if (!vertex_degree_map->count(v)) {
                                vertex_degree_map->insert({v, compute_left_core_degree(G, vertex_index_map,
                                                                                         previous_candidate_map, v,
                                                                                         k - 1,
                                                                                         h)});
                            }
                            if (vertex_degree_map->at(v) >= k) {
                                ++u_degree;
                                if (!visited_set->count(v) && !vertex_set->count(v)) {
                                    u_set->insert(v);
                                }
                            }
                        }
                    }
                }
                if (u_degree >= k) {
                    candidate_map->insert({u, u_degree});
                    next_vertex_set->merge(*u_set);
                } else {
                    remove_unsatisfied_vertices(G, u, candidate_map, evicted_set,  k, h);
                }
            }
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return true;
    }

    bool branch_multiple_core_maintenance::left_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t count = 0;
        for (auto iter = affected_set->begin(); iter != affected_set->end();) {
            auto u = *iter;
            ++iter;
            if (vertex_index_map->at(u)->count(k - 1, h) || previous_candidate_map->count(u)) {
                if (!vertex_index_map->at(u)->count(k, h)) {
                    vertex_degree_map->insert({u, compute_left_core_degree(vertex_edge_index_map,
                                                                            vertex_index_map,
                                                                            previous_candidate_map, u, k - 1, h)});

                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }else{
                        ++count;
                    }
                }
            } else {
               ++count;
            }
        }

        if (count == affected_set->size()) {
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                uint32_t u_degree = 0;
                auto u_set = make_shared<unordered_set<uint32_t>>();

                auto u_map = vertex_edge_index_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if ((vertex_index_map->at(v)->count(k - 1, h) || previous_candidate_map->count(v)) &&
                            !evicted_set->count(v)) {
                            if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v)) {
                                ++u_degree;
                            } else {
                                if (!vertex_degree_map->count(v)) {
                                    vertex_degree_map->insert({v, compute_left_core_degree(vertex_edge_index_map,
                                                                                        vertex_index_map,
                                                                                        previous_candidate_map, v,
                                                                                             k - 1, h)});
                                }

                                if (vertex_degree_map->at(v) >= k) {
                                    ++u_degree;
                                    if (!visited_set->count(v) && !vertex_set->count(v)) {
                                        u_set->insert(v);
                                    }
                                }
                            }
                        }
                    }
                }
                if (u_degree >= k) {
                    candidate_map->insert({u, u_degree});
                    next_vertex_set->merge(*u_set);
                } else {
                    remove_unsatisfied_vertices(vertex_edge_index_map, u, candidate_map, evicted_set,
                                                k, h);
                }
            }
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return true;
    }

    bool branch_multiple_core_maintenance::left_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                uint32_t k,
                                                                uint32_t h,
                                                                const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<shared_mutex>();

        {
            auto location_vector = pool->split_task(previous_candidate_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &[u, u_degree] = *iter;
                        vertex_degree_map->at(u) = u_degree;
                        sub_clear_set->insert(u);
                    }

                    global_mutex->lock();
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }


        auto count = make_shared<uint32_t>(0);
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();
                    uint32_t sub_count = 0;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->count(k - 1, h) || previous_candidate_map->count(u)) {
                            if (!vertex_index_map->at(u)->count(k, h)) {
                                if(vertex_degree_map->at(u) == 0){
                                    sub_clear_set->insert(u);
                                    vertex_degree_map->at(u) = compute_left_core_degree(G, vertex_index_map,
                                                                                        previous_candidate_map,
                                                                                        u, k - 1, h);
                                }

                                if (vertex_degree_map->at(u) >= k) {
                                    sub_vertex_set->insert(u);
                                }else{
                                    ++sub_count;
                                }
                            }
                        }else{
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

        if (*count == affected_set->size()) {
            assign(clear_set, vertex_degree_map, 0, pool);
            return false;
        }


        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto invalid_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_invalid_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        uint32_t u_degree = 0;
                        auto u_set = make_shared<unordered_set<uint32_t>>();
                        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if ((vertex_index_map->at(v)->count(k - 1, h) || previous_candidate_map->count(v)) &&
                                e_set->size() >= h && !evicted_set->count(v)) {

                                if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v) ||
                                    sub_candidate_map->count(v)) {
                                    ++u_degree;
                                } else {
                                    vertex_mutex_map->at(v)->lock();
                                    if (vertex_degree_map->at(v) == 0) {
                                        sub_clear_set->insert(v);
                                        vertex_degree_map->at(v) = compute_left_core_degree(G, vertex_index_map,
                                                                                            previous_candidate_map,
                                                                                            v, k - 1,
                                                                                            h);
                                    }
                                    vertex_mutex_map->at(v)->unlock();

                                    if (vertex_degree_map->at(v) >= k) {
                                        ++u_degree;
                                        if (!visited_set->count(v) && !vertex_set->count(v)) {
                                            u_set->insert(v);
                                        }
                                    }
                                }
                            }
                        }
                        if (u_degree >= k) {
                            sub_candidate_map->insert({u, u_degree});
                            sub_next_vertex_set->merge(*u_set);
                        } else {
                            sub_invalid_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_map->merge(*sub_candidate_map);
                    next_vertex_set->merge(*sub_next_vertex_set);
                    invalid_set->merge(*sub_invalid_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            pool->submit_task([=]{
                candidate_map->merge(*layer_candidate_map);
            });
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
            pool->barrier();
            remove_unsatisfied_vertices(G, vertex_mutex_map, invalid_set,  candidate_map, evicted_set,k, h, pool);
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }


    bool branch_multiple_core_maintenance::left_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<shared_mutex>();

        {
            auto location_vector = pool->split_task(previous_candidate_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &[u, u_degree] = *iter;
                        vertex_degree_map->at(u) = u_degree;
                        sub_clear_set->insert(u);
                    }

                    global_mutex->lock();
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto count = make_shared<uint32_t>(0);
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();
                    uint32_t sub_count = 0;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->count(k - 1, h) || previous_candidate_map->count(u)) {
                            if (!vertex_index_map->at(u)->count(k, h)) {
                                if(vertex_degree_map->at(u) == 0){
                                    sub_clear_set->insert(u);
                                    vertex_degree_map->at(u) = compute_left_core_degree(vertex_edge_index_map,
                                                                                        vertex_index_map,
                                                                                        previous_candidate_map,
                                                                                        u, k - 1, h);
                                }

                                if (vertex_degree_map->at(u) >= k) {
                                    sub_vertex_set->insert(u);
                                }else{
                                    ++sub_count;
                                }
                            }
                        } else {
                            ++sub_count;
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    *count += sub_count;
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        if (*count == affected_set->size()) {
            assign(clear_set, vertex_degree_map, 0, pool);
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto invalid_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_invalid_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                        auto &u = *iter1;
                        uint32_t u_degree = 0;
                        auto u_set = make_shared<unordered_set<uint32_t>>();
                        auto &u_map = vertex_edge_index_map->at(u);
                        for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                            for (const auto &v: *iter2->second) {
                                if ((vertex_index_map->at(v)->count(k - 1, h) || previous_candidate_map->count(v)) &&
                                    !evicted_set->count(v)) {
                                    if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v) ||
                                        sub_candidate_map->count(v)) {
                                        ++u_degree;
                                    } else {
                                        vertex_mutex_map->at(v)->lock();
                                        if (vertex_degree_map->at(v) == 0) {
                                            sub_clear_set->insert(v);
                                            vertex_degree_map->at(v) = compute_left_core_degree(
                                                    vertex_edge_index_map, vertex_index_map,
                                                    previous_candidate_map,
                                                    v, k - 1,
                                                    h);
                                        }
                                        vertex_mutex_map->at(v)->unlock();

                                        if (vertex_degree_map->at(v) >= k) {
                                            ++u_degree;
                                            if (!visited_set->count(v) && !vertex_set->count(v)) {
                                                u_set->insert(v);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (u_degree >= k) {
                            sub_candidate_map->insert({u, u_degree});
                            sub_vertex_set->merge(*u_set);
                        } else {
                            sub_invalid_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_map->merge(*sub_candidate_map);
                    next_vertex_set->merge(*sub_vertex_set);
                    invalid_set->merge(*sub_invalid_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            pool->submit_task([=]{
                candidate_map->merge(*layer_candidate_map);
            });
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
            pool->barrier();

            remove_unsatisfied_vertices(vertex_edge_index_map, vertex_mutex_map, invalid_set,
                                        candidate_map, evicted_set, k, h, pool);
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }


    bool branch_multiple_core_maintenance::middle_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                  const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                  uint32_t k,
                                                                  uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t  count = 0;
        for (auto iter = affected_set->begin(); iter != affected_set->end();) {
            auto u = *iter;
            ++iter;
            if (vertex_index_map->at(u)->count(k - 1, h - 1) || previous_candidate_map->count(u)) {
                if (!vertex_index_map->at(u)->count(k, h)) {
                    vertex_degree_map->insert({u, compute_middle_core_degree(G, vertex_index_map,
                                                                              previous_candidate_map, u,
                                                                                   k - 1, h - 1)});
                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }else{
                        ++count;
                    }
                }
            } else {
                ++count;
            }
        }

        if (count == affected_set->size()) {
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                uint32_t u_degree = 0;
                auto u_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if ((vertex_index_map->at(v)->count(k - 1, h - 1) || previous_candidate_map->count(v)) &&
                        e_set->size() >= h && !evicted_set->count(v)) {
                        if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v)) {
                            ++u_degree;
                        } else {
                            if (!vertex_degree_map->count(v)) {
                                vertex_degree_map->insert({v, compute_middle_core_degree(G, vertex_index_map,
                                                                                      previous_candidate_map, v,
                                                                                           k - 1, h - 1)});
                            }
                            if (vertex_degree_map->at(v) >= k) {
                                ++u_degree;
                                if (!visited_set->count(v) && !vertex_set->count(v)) {
                                    u_set->insert(v);
                                }
                            }
                        }
                    }
                }
                if (u_degree >= k) {
                    candidate_map->insert({u, u_degree});
                    next_vertex_set->merge(*u_set);
                } else {
                    remove_unsatisfied_vertices(G, u, candidate_map, evicted_set, k, h);
                }
            }
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return true;
    }

    bool
    branch_multiple_core_maintenance::middle_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t  count = 0;
        for (auto iter = affected_set->begin(); iter != affected_set->end();) {
            auto u = *iter;
            ++iter;
            if (vertex_index_map->at(u)->count(k - 1, h - 1) || previous_candidate_map->count(u)) {
                if (!vertex_index_map->at(u)->count(k, h)) {
                    vertex_degree_map->insert({u, compute_middle_core_degree(vertex_edge_index_map,
                                                                              vertex_index_map,
                                                                              previous_candidate_map, u, k - 1,
                                                                                   h - 1)});

                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }else{
                        ++count;
                    }
                }
            } else {
               ++count;
            }
        }

        if (count == affected_set->size()) {
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                uint32_t u_degree = 0;
                auto u_set = make_shared<unordered_set<uint32_t>>();

                auto &u_map = vertex_edge_index_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if ((vertex_index_map->at(v)->count(k - 1, h - 1) || previous_candidate_map->count(v)) &&
                            !evicted_set->count(v)) {
                            if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v)) {
                                ++u_degree;
                            } else {
                                if (!vertex_degree_map->count(v)) {
                                    vertex_degree_map->insert({v, compute_middle_core_degree(vertex_edge_index_map,
                                                                                          vertex_index_map,
                                                                                          previous_candidate_map,
                                                                                          v, k - 1, h - 1)});
                                }

                                if (vertex_degree_map->at(v) >= k) {
                                    ++u_degree;
                                    if (!visited_set->count(v) && !vertex_set->count(v)) {
                                        u_set->insert(v);
                                    }
                                }
                            }
                        }
                    }
                }

                if (u_degree >= k) {
                    candidate_map->insert({u, u_degree});
                    next_vertex_set->merge(*u_set);
                } else {
                    remove_unsatisfied_vertices(vertex_edge_index_map, u, candidate_map, evicted_set,
                                                k, h);
                }
            }
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return true;
    }

    bool branch_multiple_core_maintenance::middle_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                  const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                  uint32_t k,
                                                                  uint32_t h,
                                                                  const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<shared_mutex>();
        {
            auto location_vector = pool->split_task(previous_candidate_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &[u, u_degree] = *iter;
                        vertex_degree_map->at(u) = u_degree;
                        sub_clear_set->insert(u);
                    }

                    global_mutex->lock();
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto count = make_shared<uint32_t>(0);
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();
                    uint32_t sub_count = 0;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if ((vertex_index_map->at(u)->count(k - 1, h - 1) || previous_candidate_map->count(u))) {
                            if (!vertex_index_map->at(u)->count(k, h)) {
                                if(vertex_degree_map->at(u) == 0){
                                    sub_clear_set->insert(u);
                                    vertex_degree_map->at(u) = compute_middle_core_degree(G,
                                                                                          vertex_index_map,
                                                                                          previous_candidate_map,
                                                                                          u, k - 1, h - 1);
                                }

                                if (vertex_degree_map->at(u) >= k) {
                                    sub_vertex_set->insert(u);
                                }else{
                                    ++sub_count;
                                }
                            }
                        } else {
                            ++sub_count;
                        }
                    }
                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    *count += sub_count;
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }


        if (*count == affected_set->size()) {
            assign(clear_set, vertex_degree_map, 0, pool);
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();

        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto invalid_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_invalid_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        uint32_t u_degree = 0;
                        auto u_set = make_shared<unordered_set<uint32_t>>();
                        for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if ((vertex_index_map->at(v)->count(k - 1, h - 1) ||
                                 previous_candidate_map->count(v)) && v_edge_set->size() >= h &&
                                !evicted_set->count(v)) {
                                if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v) ||
                                    sub_candidate_map->count(v)) {
                                    ++u_degree;
                                } else {
                                    vertex_mutex_map->at(v)->lock();
                                    if (vertex_degree_map->at(v) == 0) {
                                        sub_clear_set->insert(v);
                                        vertex_degree_map->at(v) = compute_middle_core_degree(
                                                G, vertex_index_map,
                                                previous_candidate_map,
                                                v, k - 1, h - 1);
                                    }
                                    vertex_mutex_map->at(v)->unlock();

                                    if (vertex_degree_map->at(v) >= k) {
                                        ++u_degree;
                                        if (!visited_set->count(v) && !vertex_set->count(v)) {
                                            u_set->insert(v);
                                        }
                                    }
                                }
                            }
                        }
                        if (u_degree >= k) {
                            sub_candidate_map->insert({u, u_degree});
                            sub_vertex_set->merge(*u_set);
                        } else {
                            sub_invalid_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_map->merge(*sub_candidate_map);
                    next_vertex_set->merge(*sub_vertex_set);
                    invalid_set->merge(*sub_invalid_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            pool->submit_task([=]{
                candidate_map->merge(*layer_candidate_map);
            });
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
            pool->barrier();

            remove_unsatisfied_vertices(G, vertex_mutex_map, invalid_set,
                                        candidate_map, evicted_set, k, h, pool);

        }

        assign(clear_set, vertex_degree_map, 0, pool);

        return true;
    }

    bool branch_multiple_core_maintenance::middle_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<shared_mutex>();

        {
            auto location_vector = pool->split_task(previous_candidate_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &[u, u_degree] = *iter;
                        vertex_degree_map->at(u) = u_degree;
                        sub_clear_set->insert(u);
                    }

                    global_mutex->lock();
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto count = make_shared<uint32_t>(0);
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();
                    uint32_t sub_count = 0;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if ((vertex_index_map->at(u)->count(k - 1, h - 1) || previous_candidate_map->count(u))) {
                            if (!vertex_index_map->at(u)->count(k, h)) {
                                if(vertex_degree_map->at(u) == 0){
                                    sub_clear_set->insert(u);
                                    vertex_degree_map->at(u) = compute_middle_core_degree(vertex_edge_index_map,
                                                                                          vertex_index_map,
                                                                                          previous_candidate_map,
                                                                                          u, k - 1, h - 1);
                                }

                                if (vertex_degree_map->at(u) >= k) {
                                    sub_vertex_set->insert(u);
                                }else{
                                    ++sub_count;
                                }
                            }
                        } else {
                            ++sub_count;
                        }
                    }
                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    *count += sub_count;
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }


        if (*count == affected_set->size()) {
            assign(clear_set, vertex_degree_map, 0, pool);
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();

        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto invalid_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_invalid_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        uint32_t u_degree = 0;
                        auto u_set = make_shared<unordered_set<uint32_t>>();
                        auto &u_map = vertex_edge_index_map->at(u);
                        for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                            for (const auto &v: *iter2->second) {
                                if ((vertex_index_map->at(v)->count(k - 1, h - 1) ||
                                     previous_candidate_map->count(v)) && !evicted_set->count(v)) {
                                    if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v) ||
                                        sub_candidate_map->count(v)) {
                                        ++u_degree;
                                    } else {
                                        vertex_mutex_map->at(v)->lock();
                                        if (vertex_degree_map->at(v) == 0) {
                                            sub_clear_set->insert(v);
                                            vertex_degree_map->at(v) = compute_middle_core_degree(
                                                    vertex_edge_index_map, vertex_index_map,
                                                    previous_candidate_map,
                                                    v, k - 1, h - 1);
                                        }
                                        vertex_mutex_map->at(v)->unlock();

                                        if (vertex_degree_map->at(v) >= k) {
                                            ++u_degree;
                                            if (!visited_set->count(v) && !vertex_set->count(v)) {
                                                u_set->insert(v);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (u_degree >= k) {
                            sub_candidate_map->insert({u, u_degree});
                            sub_vertex_set->merge(*u_set);
                        } else {
                            sub_invalid_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_map->merge(*sub_candidate_map);
                    next_vertex_set->merge(*sub_vertex_set);
                    invalid_set->merge(*sub_invalid_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            pool->submit_task([=]{
                candidate_map->merge(*layer_candidate_map);
            });
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
            pool->barrier();

            remove_unsatisfied_vertices(vertex_edge_index_map, vertex_mutex_map, invalid_set,
                                        candidate_map, evicted_set, k, h, pool);
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }


    bool branch_multiple_core_maintenance::right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                 const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                 uint32_t k,
                                                                 uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t  count = 0;
        for (const auto &u:*affected_set) {
            if (vertex_index_map->at(u)->count(k, h - 1) || previous_candidate_map->count(u)) {
                if (!vertex_index_map->at(u)->count(k, h)) {
                    vertex_degree_map->insert({u, compute_right_core_degree(G, vertex_index_map,
                                                                             previous_candidate_map, u, k, h - 1)});

                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }else{
                        ++count;
                    }
                }
            }else{
                ++count;
            }
        }

        if (count == affected_set->size()) {
            return false;
        }

        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                uint32_t u_degree = 0;
                auto u_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if ((vertex_index_map->at(v)->count(k, h - 1) || previous_candidate_map->count(v)) &&
                        e_set->size() >= h && !evicted_set->count(v)) {
                        if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v)) {
                            ++u_degree;
                        } else {
                            if (!vertex_degree_map->count(v)) {
                                vertex_degree_map->insert({v, compute_right_core_degree(G, vertex_index_map,
                                                                                     previous_candidate_map, v, k,
                                                                                          h - 1)});
                            }
                            if (vertex_degree_map->at(v) >= k) {
                                ++u_degree;
                                if (!visited_set->count(v) && !vertex_set->count(v)) {
                                    u_set->insert(v);
                                }
                            }
                        }
                    }
                }
                if (u_degree >= k) {
                    candidate_map->insert({u, u_degree});
                    next_vertex_set->merge(*u_set);
                } else {
                    remove_unsatisfied_vertices(G, u, candidate_map, evicted_set, k, h);
                }
            }
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return true;
    }

    bool
    branch_multiple_core_maintenance::right_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        uint32_t count = 0;
        for (const auto &u:*affected_set) {
            if (vertex_index_map->at(u)->count(k, h - 1) || previous_candidate_map->count(u)) {
                if (!vertex_index_map->at(u)->count(k, h)) {
                    vertex_degree_map->insert({u, compute_right_core_degree(vertex_edge_index_map,
                                                                             vertex_index_map,
                                                                             previous_candidate_map, u, k, h - 1)});

                    if (vertex_degree_map->at(u) >= k) {
                        vertex_set->insert(u);
                    }else{
                        ++count;
                    }
                }
            } else {
                ++count;
            }
        }

        if (count == affected_set->size()) {
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                uint32_t u_degree = 0;
                auto u_set = make_shared<unordered_set<uint32_t>>();

                auto &u_map = vertex_edge_index_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if ((vertex_index_map->at(v)->count(k, h - 1) || previous_candidate_map->count(v)) &&
                            !evicted_set->count(v)) {
                            if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v)) {
                                ++u_degree;
                            } else {
                                if (!vertex_degree_map->count(v)) {
                                    vertex_degree_map->insert({v, compute_right_core_degree(vertex_edge_index_map,
                                                                                         vertex_index_map,
                                                                                         previous_candidate_map, v,
                                                                                         k, h - 1)});
                                }

                                if (vertex_degree_map->at(v) >= k) {
                                    ++u_degree;
                                    if (!visited_set->count(v) && !vertex_set->count(v)) {
                                        u_set->insert(v);
                                    }
                                }
                            }
                        }
                    }
                }

                if (u_degree >= k) {
                    candidate_map->insert({u, u_degree});
                    next_vertex_set->merge(*u_set);
                } else {
                    remove_unsatisfied_vertices(vertex_edge_index_map, u, candidate_map, evicted_set, 
                                                k, h);
                }
            }
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return true;
    }

    bool branch_multiple_core_maintenance::right_candidate_graph(const shared_ptr<temporal_graph> &G,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                 const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                 uint32_t k,
                                                                 uint32_t h,
                                                                 const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<shared_mutex>();

        {
            auto location_vector = pool->split_task(previous_candidate_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &[u, u_degree] = *iter;
                        vertex_degree_map->at(u) = u_degree;
                        sub_clear_set->insert(u);
                    }

                    global_mutex->lock();
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto count = make_shared<uint32_t>(0);
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();
                    uint32_t sub_count = 0;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if ((vertex_index_map->at(u)->count(k, h - 1) || previous_candidate_map->count(u))) {
                            if (!vertex_index_map->at(u)->count(k, h)) {
                                if(vertex_degree_map->at(u) == 0){
                                    sub_clear_set->insert(u);
                                    vertex_degree_map->at(u) = compute_right_core_degree(G, vertex_index_map,
                                                                                         previous_candidate_map,
                                                                                         u, k, h - 1);
                                }

                                if (vertex_degree_map->at(u) >= k) {
                                    sub_vertex_set->insert(u);
                                }else{
                                    ++sub_count;
                                }
                            }
                        } else {
                            ++sub_count;
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    *count += sub_count;
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }


        if (*count == affected_set->size()) {
            assign(clear_set, vertex_degree_map, 0, pool);
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto invalid_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_invalid_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        uint32_t u_degree = 0;
                        auto u_set = make_shared<unordered_set<uint32_t>>();
                        for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if ((vertex_index_map->at(v)->count(k, h - 1) || previous_candidate_map->count(v)) &&
                                e_set->size() >= h && !evicted_set->count(v)) {
                                if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v) ||
                                    sub_candidate_map->count(v)) {
                                    ++u_degree;
                                } else {
                                    vertex_mutex_map->at(v)->lock();
                                    if (vertex_degree_map->at(v) == 0) {
                                        sub_clear_set->insert(v);
                                        vertex_degree_map->at(v) = compute_right_core_degree(G, vertex_index_map,
                                                                                             previous_candidate_map,
                                                                                             v, k,
                                                                                                  h - 1);
                                    }
                                    vertex_mutex_map->at(v)->unlock();

                                    if (vertex_degree_map->at(v) >= k) {
                                        ++u_degree;
                                        if (!visited_set->count(v) && !vertex_set->count(v)) {
                                            u_set->insert(v);
                                        }
                                    }
                                }
                            }
                        }
                        if (u_degree >= k) {
                            sub_candidate_map->insert({u, u_degree});
                            sub_vertex_set->merge(*u_set);
                        } else {
                            sub_invalid_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_map->merge(*sub_candidate_map);
                    next_vertex_set->merge(*sub_vertex_set);
                    invalid_set->merge(*sub_invalid_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();

            pool->submit_task([=]{
                candidate_map->merge(*layer_candidate_map);
            });
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
            pool->barrier();
            remove_unsatisfied_vertices(G, vertex_mutex_map, invalid_set, candidate_map, evicted_set,
                                        k, h, pool);
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }

    bool branch_multiple_core_maintenance::right_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_candidate_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<shared_mutex>();

        {
            auto location_vector = pool->split_task(previous_candidate_map);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        const auto &[u, u_degree] = *iter;
                        vertex_degree_map->at(u) = u_degree;
                        sub_clear_set->insert(u);
                    }

                    global_mutex->lock();
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto count = make_shared<uint32_t>(0);
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();
                    uint32_t sub_count = 0;

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->count(k, h - 1) || previous_candidate_map->count(u)) {
                            if (!vertex_index_map->at(u)->count(k, h)) {
                                if(vertex_degree_map->at(u) == 0){
                                    sub_clear_set->insert(u);
                                    vertex_degree_map->at(u) = compute_right_core_degree(vertex_edge_index_map,
                                                                                         vertex_index_map,
                                                                                         previous_candidate_map,
                                                                                         u, k, h - 1);
                                }

                                if (vertex_degree_map->at(u) >= k) {
                                    sub_vertex_set->insert(u);
                                } else {
                                    ++sub_count;
                                }
                            }
                        } else {
                            ++sub_count;
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    *count += sub_count;
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        if (*count == affected_set->size()) {
            assign(clear_set, vertex_degree_map, 0, pool);
            return false;
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto invalid_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            auto layer_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto location_vector = pool->split_task(vertex_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_invalid_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                        auto &u = *iter1;
                        uint32_t u_degree = 0;
                        auto u_set = make_shared<unordered_set<uint32_t>>();
                        auto &u_map = vertex_edge_index_map->at(u);
                        for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                            for (const auto &v: *iter2->second) {
                                if ((vertex_index_map->at(v)->count(k, h - 1) || previous_candidate_map->count(v)) &&
                                    !evicted_set->count(v)) {
                                    if (vertex_index_map->at(v)->count(k, h) || candidate_map->count(v) ||
                                        sub_candidate_map->count(v)) {
                                        ++u_degree;
                                    } else {
                                        vertex_mutex_map->at(v)->lock();
                                        if (vertex_degree_map->at(v) == 0) {
                                            sub_clear_set->insert(v);
                                            vertex_degree_map->at(v) = compute_right_core_degree(
                                                    vertex_edge_index_map, vertex_index_map,
                                                    previous_candidate_map,
                                                    v, k,
                                                    h - 1);
                                        }
                                        vertex_mutex_map->at(v)->unlock();

                                        if (vertex_degree_map->at(v) >= k) {
                                            ++u_degree;
                                            if (!visited_set->count(v) && !vertex_set->count(v)) {
                                                u_set->insert(v);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (u_degree >= k) {
                            sub_candidate_map->insert({u, u_degree});
                            sub_vertex_set->merge(*u_set);
                        } else {
                            sub_invalid_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    layer_candidate_map->merge(*sub_candidate_map);
                    next_vertex_set->merge(*sub_vertex_set);
                    invalid_set->merge(*sub_invalid_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
            pool->submit_task([=]{
                candidate_map->merge(*layer_candidate_map);
            });
            visited_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
            pool->barrier();

            remove_unsatisfied_vertices(vertex_edge_index_map, vertex_mutex_map, invalid_set,
                                        candidate_map, evicted_set, k, h, pool);
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return true;
    }


    void branch_multiple_core_maintenance::left_decomposition(const shared_ptr<temporal_graph> &G,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                              uint32_t k,
                                                              uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        while (!candidate_map->empty()) {
            for (const auto &[u, u_degree]: *candidate_map) {
                if (u_degree < k) {
                    vertex_set->insert(u);
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    vertex_index_map->at(u)->insert(k - 1, h);
                }
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                    vertex_index_map->at(u)->insert(k - 1, h);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
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
            ++k;
        }
    }

    void branch_multiple_core_maintenance::left_decomposition(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t> >();
        while (!candidate_map->empty()) {
            for (const auto &[u, u_degree]: *candidate_map) {
                if (u_degree < k) {
                    vertex_set->insert(u);
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    vertex_index_map->at(u)->insert(k - 1, h);
                }
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                    vertex_index_map->at(u)->insert(k - 1, h);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
                    auto edge_map = vertex_edge_index_map->at(u);
                    for (auto iter = edge_map->lower_bound(h); iter != edge_map->end(); ++iter) {
                        for (const auto &v: *iter->second) {
                            if (candidate_map->count(v) && candidate_map->at(v) >= k) {
                                --candidate_map->at(v);
                                if (candidate_map->at(v) < k) {
                                    next_vertex_set->insert(v);
                                }
                            }
                        }
                    }
                }
                swap(*vertex_set, *next_vertex_set);
            }
            ++k;
        }
    }

    void branch_multiple_core_maintenance::left_decomposition(const shared_ptr<temporal_graph> &G,
                                                              const shared_ptr<mutex> &global_mutex,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                              uint32_t k,
                                                              uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t> >();
        while (!candidate_map->empty()) {
            for (const auto &[u, u_degree]: *candidate_map) {
                if (u_degree < k) {
                    vertex_set->insert(u);
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    new_vertex_index_map->at(u)->insert(k - 1, h);
                }
                global_mutex->unlock();
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(k - 1, h);
                }
                global_mutex->unlock();
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
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
            ++k;
        }
    }

    void branch_multiple_core_maintenance::left_decomposition(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<mutex> &global_mutex,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t> >();
        while (!candidate_map->empty()) {
            for (const auto &[u, u_degree]: *candidate_map) {
                if (u_degree < k) {
                    vertex_set->insert(u);
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    new_vertex_index_map->at(u)->insert(k - 1, h);
                }
                global_mutex->unlock();
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(k - 1, h);
                }
                global_mutex->unlock();
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
                    auto edge_map = vertex_edge_index_map->at(u);
                    for (auto iter = edge_map->lower_bound(h); iter != edge_map->end(); ++iter) {
                        for (const auto &v: *iter->second) {
                            if (candidate_map->count(v) && candidate_map->at(v) >= k) {
                                --candidate_map->at(v);
                                if (candidate_map->at(v) < k) {
                                    next_vertex_set->insert(v);
                                }
                            }
                        }
                    }
                }
                swap(*vertex_set, *next_vertex_set);
            }
            ++k;
        }
    }

    void branch_multiple_core_maintenance::left_decomposition(const shared_ptr<temporal_graph> &G,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                              uint32_t k,
                                                              uint32_t h,
                                                              const shared_ptr<thread_pool> &pool) {
        auto global_mutex = make_shared<mutex>();
        auto thread_number = pool->get_thread_number();

        auto vertex_set = make_shared<unordered_set<uint32_t> >();
        while (!candidate_map->empty()) {
            {
                auto location_vector = pool->split_task(candidate_map);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &[u, u_degree] = *iter;
                            if (u_degree < k) {
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


            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                {
                    auto location_vector = pool->split_task(vertex_set);
                    for(uint32_t i = 0; i < thread_number; ++i){
                        pool->submit_task([=]{
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &u = *iter;
                                if(!new_vertex_index_map->at(u)){
                                    new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                }
                                new_vertex_index_map->at(u)->insert(k - 1, h);
                            }
                        });
                    }
                    candidate_map->clear();
                    pool->barrier();
                }
                break;
            }

            {
                while (!vertex_set->empty()) {
                    for (const auto &u: *vertex_set) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }
                    {
                        auto location_vector = pool->split_task(vertex_set);
                        for(uint32_t i = 0; i < thread_number; ++i){
                            pool->submit_task([=]{
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                    auto &u = *iter;
                                    if(!new_vertex_index_map->at(u)){
                                        new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                    }
                                    new_vertex_index_map->at(u)->insert(k - 1, h);
                                }
                            });
                        }
                        for(const auto &u: *vertex_set){
                            candidate_map->erase(u);
                        }
                        pool->barrier();
                    }
                    auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto location_vector = pool->split_task(vertex_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();

                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                auto &u = *iter;
                                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                    if (candidate_map->count(v) && e_set->size() >= h  && candidate_map->at(v) >= k) {
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
            ++k;
        }
    }

    void branch_multiple_core_maintenance::left_decomposition(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t> >();
        auto global_mutex = make_shared<mutex>();

        auto thread_number = pool->get_thread_number();
        while (!candidate_map->empty()) {
            {
                auto location_vector = pool->split_task(candidate_map);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &[u, u_degree] = *iter;
                            if (u_degree < k) {
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

            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                {
                    auto location_vector = pool->split_task(vertex_set);
                    for(uint32_t i = 0; i < thread_number; ++i){
                        pool->submit_task([=]{
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &u = *iter;
                                if(!new_vertex_index_map->at(u)){
                                    new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                }
                                new_vertex_index_map->at(u)->insert(k - 1, h);
                            }
                        });
                    }
                    candidate_map->clear();
                    pool->barrier();
                }
                break;
            }

            {
                while (!vertex_set->empty()) {
                    for (const auto &u: *vertex_set) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }
                    {
                        auto location_vector = pool->split_task(vertex_set);
                        for(uint32_t i = 0; i < thread_number; ++i){
                            pool->submit_task([=]{
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                    auto &u = *iter;
                                    if(!new_vertex_index_map->at(u)){
                                        new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                    }
                                    new_vertex_index_map->at(u)->insert(k - 1, h);
                                }
                            });
                        }
                        for(const auto &u: *vertex_set){
                            candidate_map->erase(u);
                        }
                        pool->barrier();
                    }
                    auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto location_vector = pool->split_task(vertex_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                            for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                                auto &u = *iter1;
                                auto &u_map = vertex_edge_index_map->at(u);
                                for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                                    for (const auto &v: *iter2->second) {
                                        if (candidate_map->count(v) && candidate_map->at(v) >= k) {

                                            vertex_mutex_map->at(v)->lock();
                                            --candidate_map->at(v);
                                            vertex_mutex_map->at(v)->unlock();

                                            if (candidate_map->at(v) < k) {
                                                sub_next_vertex_set->insert(v);
                                            }
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
            ++k;
        }
    }

    void branch_multiple_core_maintenance::middle_decomposition(const shared_ptr<temporal_graph> &G,
                                                                const shared_ptr<mutex> &global_mutex,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                                uint32_t delta,
                                                                const shared_ptr<thread_pool> &pool) {
        while (!candidate_map->empty()) {
            {
                auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(*candidate_map);
                pool->submit_task([=] {
                    left_decomposition(G, global_mutex, sub_candidate_map, new_vertex_index_map,
                                       delta + 1, delta);
                });

            }

            {
                auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(*candidate_map);
                pool->submit_task([=] {
                    right_decomposition(G, global_mutex, sub_candidate_map, new_vertex_index_map,
                                        delta, delta + 1);
                });
            }

            auto vertex_set = make_shared<unordered_set<uint32_t>>();

            for (const auto &[u, u_degree]: *candidate_map) {
                if (candidate_map->at(u) < delta + 1) {
                    vertex_set->insert(u);
                } else {
                    for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                        if (candidate_map->count(v) && e_set->size() == delta) {
                            --candidate_map->at(u);
                        }
                    }
                    if (candidate_map->at(u) < delta + 1) {
                        vertex_set->insert(u);
                    }
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(delta, delta);
                }
                global_mutex->unlock();
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(delta, delta);
                }
                global_mutex->unlock();
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
                    for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                        if (candidate_map->count(v)  && e_set->size() >= delta + 1 && candidate_map->at(v) >= delta + 1) {
                            --candidate_map->at(v);
                            if (candidate_map->at(v) < delta + 1) {
                                next_vertex_set->insert(v);
                            }
                        }
                    }
                }
                swap(*vertex_set, *next_vertex_set);
            }
            ++delta;
        }
        pool->barrier();
    }

    void branch_multiple_core_maintenance::middle_decomposition(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<mutex> &global_mutex,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
            uint32_t delta,
            const shared_ptr<thread_pool> &pool) {
        while (!candidate_map->empty()) {
            {
                auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(*candidate_map);
                pool->submit_task([=] {
                    left_decomposition(vertex_edge_index_map, global_mutex, sub_candidate_map,
                                       new_vertex_index_map,
                                       delta + 1, delta);
                });
            }

            {
                auto sub_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(*candidate_map);
                pool->submit_task([=] {
                    right_decomposition(vertex_edge_index_map, global_mutex, sub_candidate_map,
                                        new_vertex_index_map,
                                        delta, delta + 1);
                });
            }

            auto vertex_set = make_shared<unordered_set<uint32_t>>();

            for (const auto &[u, u_degree]: *candidate_map) {
                if (u_degree < delta + 1) {
                    vertex_set->insert(u);
                } else {
                    auto edge_map = vertex_edge_index_map->at(u);
                    if (edge_map->count(delta)) {
                        for (const auto &v: *edge_map->at(delta)) {
                            if (candidate_map->count(v)) {
                                --candidate_map->at(u);
                            }
                        }
                    }
                    if (candidate_map->at(u) < delta + 1) {
                        vertex_set->insert(u);
                    }
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(delta, delta);
                }
                global_mutex->unlock();
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(delta, delta);
                }
                global_mutex->unlock();
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
                    auto &edge_map = vertex_edge_index_map->at(u);
                    for (auto iter = edge_map->lower_bound(delta + 1); iter != edge_map->end(); ++iter) {
                        for (const auto &v: *iter->second) {
                            if (candidate_map->count(v) && candidate_map->at(v) >= delta + 1) {
                                --candidate_map->at(v);
                                if (candidate_map->at(v) < delta + 1) {
                                    next_vertex_set->insert(v);
                                }
                            }
                        }
                    }
                }
                swap(*vertex_set, *next_vertex_set);
            }
            ++delta;
        }
        pool->barrier();
    }

    void branch_multiple_core_maintenance::middle_decomposition(const shared_ptr<temporal_graph> &G,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                                uint32_t delta,
                                                                const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        while (!candidate_map->empty()) {
            auto sub_candidate_map1 = make_shared<unordered_map<uint32_t, uint32_t>>();
            pool->submit_task([=]{
               copy(candidate_map->begin(), candidate_map->end(), std::inserter(*sub_candidate_map1, sub_candidate_map1->end()));
            });
            auto sub_candidate_map2 = make_shared<unordered_map<uint32_t, uint32_t>>(*candidate_map);
            pool->barrier();

            {
                auto &sub_candidate_map = sub_candidate_map1;
                left_decomposition(G, vertex_mutex_map, sub_candidate_map, new_vertex_index_map, delta + 1,
                                   delta, pool);
            }

            {
                auto &sub_candidate_map = sub_candidate_map2;
                right_decomposition(G, vertex_mutex_map, sub_candidate_map, new_vertex_index_map, delta,
                                    delta + 1, pool);
            }

            auto vertex_set = make_shared<unordered_set<uint32_t>>();
            {
                auto location_vector = pool->split_task(candidate_map);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &[u, u_degree] = *iter;
                            if (u_degree < delta + 1) {
                                sub_vertex_set->insert(u);
                            } else {
                                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                    if (candidate_map->count(v) && e_set->size() == delta) {
                                        --candidate_map->at(u);
                                    }
                                }
                                if (candidate_map->at(u) < delta + 1) {
                                    sub_vertex_set->insert(u);
                                }
                            }
                        }

                        global_mutex->lock();
                        vertex_set->merge(*sub_vertex_set);
                        global_mutex->unlock();
                    });
                    pool->barrier();
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                {
                    auto location_vector = pool->split_task(vertex_set);
                    for(uint32_t i = 0; i < thread_number; ++i){
                        pool->submit_task([=]{
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &u = *iter;
                                if(!new_vertex_index_map->at(u)){
                                    new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                }
                                new_vertex_index_map->at(u)->insert(delta, delta);
                            }
                        });
                    }
                    candidate_map->clear();
                    pool->barrier();
                }
                break;
            }

            {
                while (!vertex_set->empty()) {
                    for (const auto &u: *vertex_set) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }
                    {
                        auto location_vector = pool->split_task(vertex_set);
                        for(uint32_t i = 0; i < thread_number; ++i){
                            pool->submit_task([=]{
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                    auto &u = *iter;
                                    if(!new_vertex_index_map->at(u)){
                                        new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                    }
                                    new_vertex_index_map->at(u)->insert(delta, delta);
                                }
                            });
                        }
                        for(const auto &u:*vertex_set){
                            candidate_map->erase(u);
                        }
                        pool->barrier();
                    }
                    auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto location_vector = pool->split_task(vertex_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                auto &u = *iter;
                                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                    if (candidate_map->count(v) && e_set->size() >= delta + 1 &&  candidate_map->at(v) >= delta + 1) {

                                        vertex_mutex_map->at(v)->lock();
                                        --candidate_map->at(v);
                                        vertex_mutex_map->at(v)->unlock();

                                        if (candidate_map->at(v) < delta + 1) {
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
            ++delta;
        }
    }

    void branch_multiple_core_maintenance::middle_decomposition(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
            uint32_t delta,
            const shared_ptr<thread_pool> &pool) {

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        while (!candidate_map->empty()) {
            auto sub_candidate_map1 = make_shared<unordered_map<uint32_t, uint32_t>>();
            pool->submit_task([=]{
                copy(candidate_map->begin(), candidate_map->end(), std::inserter(*sub_candidate_map1, sub_candidate_map1->end()));
            });
            auto sub_candidate_map2 = make_shared<unordered_map<uint32_t, uint32_t>>(*candidate_map);
            pool->barrier();
            {
                auto &sub_candidate_map = sub_candidate_map1;
                left_decomposition(vertex_edge_index_map, vertex_mutex_map, sub_candidate_map,
                                   new_vertex_index_map,
                                   delta + 1, delta, pool);
            }

            {
                auto &sub_candidate_map = sub_candidate_map2;
                right_decomposition(vertex_edge_index_map, vertex_mutex_map, sub_candidate_map,
                                    new_vertex_index_map,
                                    delta, delta + 1, pool);
            }


            auto vertex_set = make_shared<unordered_set<uint32_t>>();
            {
                auto location_vector = pool->split_task(candidate_map);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &[u, u_degree] = *iter;
                            if (u_degree < delta + 1) {
                                sub_vertex_set->insert(u);
                            } else {
                                auto u_map = vertex_edge_index_map->at(u);
                                if (u_map->count(delta)) {
                                    for (const auto &v: *u_map->at(delta)) {
                                        if (candidate_map->count(v)) {
                                            --candidate_map->at(u);
                                        }
                                    }
                                }
                                if (candidate_map->at(u) < delta + 1) {
                                    sub_vertex_set->insert(u);
                                }
                            }
                        }

                        global_mutex->lock();
                        vertex_set->merge(*sub_vertex_set);
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
            }


            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                {
                    auto location_vector = pool->split_task(vertex_set);
                    for(uint32_t i = 0; i < thread_number; ++i){
                        pool->submit_task([=]{
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &u = *iter;
                                if(!new_vertex_index_map->at(u)){
                                    new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                }
                                new_vertex_index_map->at(u)->insert(delta, delta);
                            }
                        });
                    }
                    candidate_map->clear();
                    pool->barrier();
                }
                break;
            }

            {
                while (!vertex_set->empty()) {
                    for (const auto &u: *vertex_set) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }
                    {
                        auto location_vector = pool->split_task(vertex_set);
                        for(uint32_t i = 0; i < thread_number; ++i){
                            pool->submit_task([=]{
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                    auto &u = *iter;
                                    if(!new_vertex_index_map->at(u)){
                                        new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                    }

                                    new_vertex_index_map->at(u)->insert(delta, delta);
                                }
                            });
                        }
                        for(const auto &u:*vertex_set){
                            candidate_map->erase(u);
                        }
                        pool->barrier();
                    }
                    auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto location_vector = pool->split_task(vertex_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();

                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                                auto &u = *iter1;
                                auto &u_map = vertex_edge_index_map->at(u);
                                for (auto iter2 = u_map->lower_bound(delta + 1); iter2 != u_map->end(); ++iter2) {
                                    for (const auto &v: *iter2->second) {
                                        if (candidate_map->count(v) && candidate_map->at(v) >= delta + 1) {

                                            vertex_mutex_map->at(v)->lock();
                                            --candidate_map->at(v);
                                            vertex_mutex_map->at(v)->unlock();

                                            if (candidate_map->at(v) < delta + 1) {
                                                sub_next_vertex_set->insert(v);
                                            }
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
            ++delta;
        }
    }


    void branch_multiple_core_maintenance::right_decomposition(const shared_ptr<temporal_graph> &G,
                                                               const shared_ptr<mutex> &global_mutex,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                               uint32_t k,
                                                               uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t> >();

        while (!candidate_map->empty()) {
            for (const auto &[u, u_degree]: *candidate_map) {
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (candidate_map->count(v) && e_set->size() == h - 1) {
                        --candidate_map->at(u);
                    }
                }
                if (candidate_map->at(u) < k) {
                    vertex_set->insert(u);
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(k, h - 1);
                }
                global_mutex->unlock();
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(k, h - 1);
                }
                global_mutex->unlock();
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
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
            ++h;
        }
    }

    void branch_multiple_core_maintenance::right_decomposition(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<mutex> &global_mutex,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t> >();

        while (!candidate_map->empty()) {
            for (const auto &[u, u_degree]: *candidate_map) {
                auto &u_map = vertex_edge_index_map->at(u);
                if (u_map->count(h - 1)) {
                    for (const auto &v: *u_map->at(h - 1)) {
                        if (candidate_map->count(v)) {
                            --candidate_map->at(u);
                        }
                    }
                    if (candidate_map->at(u) < k) {
                        vertex_set->insert(u);
                    }
                }
            }

            if (vertex_set->size() == candidate_map->size()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(k, h - 1);
                }
                global_mutex->unlock();
                candidate_map->clear();
                break;
            }

            while (!vertex_set->empty()) {
                global_mutex->lock();
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    }
                    new_vertex_index_map->at(u)->insert(k, h - 1);
                }
                global_mutex->unlock();
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (const auto &u: *vertex_set) {
                    auto &u_map = vertex_edge_index_map->at(u);
                    for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                        for (const auto &v: *iter->second) {
                            if (candidate_map->count(v) && candidate_map->at(v) >= k) {
                                --candidate_map->at(v);
                                if (candidate_map->at(v) < k) {
                                    next_vertex_set->insert(v);
                                }
                            }
                        }
                    }
                }
                swap(*vertex_set, *next_vertex_set);
            }
            ++h;
        }
    }

    void branch_multiple_core_maintenance::right_decomposition(const shared_ptr<temporal_graph> &G,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
                                                               uint32_t k,
                                                               uint32_t h,
                                                               const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        while (!candidate_map->empty()) {
            auto vertex_set = make_shared<unordered_set<uint32_t> >();
            {
                auto location_vector = pool->split_task(candidate_map);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &[u, u_degree] = *iter;

                            for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                if (candidate_map->count(v) && e_set->size() == h - 1) {
                                    --candidate_map->at(u);
                                }
                            }
                            if (candidate_map->at(u) < k) {
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


            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                {
                    auto location_vector = pool->split_task(vertex_set);
                    for(uint32_t i = 0; i < thread_number; ++i){
                        pool->submit_task([=]{
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &u = *iter;
                                if(!new_vertex_index_map->at(u)){
                                    new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                }
                                new_vertex_index_map->at(u)->insert(k, h - 1);
                            }
                        });
                    }
                    candidate_map->clear();
                    pool->barrier();
                }
                break;
            }

            {
                while (!vertex_set->empty()) {
                    for (const auto &u: *vertex_set) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }
                    {
                        auto location_vector = pool->split_task(vertex_set);
                        for(uint32_t i = 0; i < thread_number; ++i){
                            pool->submit_task([=]{
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                    auto &u = *iter;
                                    if(!new_vertex_index_map->at(u)){
                                        new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                    }
                                    new_vertex_index_map->at(u)->insert(k, h - 1);
                                }
                            });
                        }
                        for(const auto &u: *vertex_set){
                            candidate_map->erase(u);
                        }
                        pool->barrier();
                    }
                    auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto location_vector = pool->split_task(vertex_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();

                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
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
            ++h;
        }
    }

    void branch_multiple_core_maintenance::right_decomposition(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &new_vertex_index_map,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        while (!candidate_map->empty()) {
            auto vertex_set = make_shared<unordered_set<uint32_t> >();
            {
                auto location_vector = pool->split_task(candidate_map);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &[u, u_degree] = *iter;
                            auto &u_map = vertex_edge_index_map->at(u);
                            if (u_map->count(h - 1)) {
                                for (const auto &v: *u_map->at(h - 1)) {
                                    if (candidate_map->count(v)) {
                                        --candidate_map->at(u);
                                    }
                                }
                                if (candidate_map->at(u) < k) {
                                    sub_vertex_set->insert(u);
                                }
                            }
                        }

                        global_mutex->lock();
                        vertex_set->merge(*sub_vertex_set);
                        global_mutex->unlock();
                    });
                }
                pool->barrier();
            }


            if (vertex_set->size() == candidate_map->size()) {
                for (const auto &u: *vertex_set) {
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                {
                    auto location_vector = pool->split_task(vertex_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                auto &u = *iter;
                                if(!new_vertex_index_map->at(u)){
                                    new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                }
                                new_vertex_index_map->at(u)->insert(k, h - 1);
                            }
                        });
                    }
                    candidate_map->clear();
                    pool->barrier();
                }
                break;
            }

            {

                while (!vertex_set->empty()) {
                    for (const auto &u: *vertex_set) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }
                    {
                        auto location_vector = pool->split_task(vertex_set);
                        for (uint32_t i = 0; i < thread_number; ++i) {
                            pool->submit_task([=] {
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                    auto &u = *iter;
                                    if(!new_vertex_index_map->at(u)){
                                        new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                    }
                                    new_vertex_index_map->at(u)->insert(k, h - 1);
                                }
                            });
                        }
                        for (const auto &u: *vertex_set) {
                            candidate_map->erase(u);
                        }
                        pool->barrier();
                    }
                    auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto location_vector = pool->split_task(vertex_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();

                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                auto &u = *iter;
                                auto &u_map = vertex_edge_index_map->at(u);
                                for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                                    for (const auto &v: *iter2->second) {
                                        if (candidate_map->count(v) && candidate_map->at(v) >= k) {

                                            vertex_mutex_map->at(v)->lock();
                                            --candidate_map->at(v);
                                            vertex_mutex_map->at(v)->unlock();

                                            if (candidate_map->at(v) < k) {
                                                sub_next_vertex_set->insert(v);
                                            }
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
            ++h;
        }
    }

    void branch_multiple_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                                       uint32_t w,
                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                                       uint32_t k,
                                                                       uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(w);

        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                evicted_set->insert(u);
                candidate_map->erase(u);
            }
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*vertex_set) {
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (candidate_map->count(v) && candidate_map->at(v) >= k && e_set->size() >= h) {
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

    void branch_multiple_core_maintenance::remove_unsatisfied_vertices(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            uint32_t w,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_set<uint32_t>> &evicted_set,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(w);

        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                evicted_set->insert(u);
                candidate_map->erase(u);
            }
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*vertex_set){
                auto &u_map = vertex_edge_index_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if (candidate_map->count(v) && candidate_map->at(v) >= k) {
                            --candidate_map->at(v);
                            if (candidate_map->at(v) < k) {
                                next_vertex_set->insert(v);
                            }
                        }
                    }
                }
            }
            swap(*vertex_set, *next_vertex_set);
        }
    }

    void branch_multiple_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                                       const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                                       uint32_t k,
                                                                       uint32_t h) {
        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                evicted_set->insert(u);
                candidate_map->erase(u);
            }
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*vertex_set) {
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (candidate_map->count(v) && candidate_map->at(v) >= k && e_set->size() >= h) {
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

    void branch_multiple_core_maintenance::remove_unsatisfied_vertices(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_set<uint32_t>> &vertex_set,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_set<uint32_t>> &evicted_set,
            uint32_t k,
            uint32_t h) {
        while (!vertex_set->empty()) {
            for(const auto &u:*vertex_set){
                evicted_set->insert(u);
                candidate_map->erase(u);
            }
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*vertex_set){
                auto &u_map = vertex_edge_index_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if (candidate_map->count(v) && candidate_map->at(v) >= k) {
                            --candidate_map->at(v);
                            if (candidate_map->at(v) < k) {
                                next_vertex_set->insert(v);
                            }
                        }
                    }
                }
            }
            swap(*vertex_set, *next_vertex_set);
        }
    }

    void branch_multiple_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<temporal_graph> &G,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &evicted_set,
                                                                       uint32_t k,
                                                                       uint32_t h,
                                                                       const shared_ptr<thread_pool> &pool) {
        auto global_mutex = make_shared<mutex>();
        auto thread_number = pool->get_thread_number();
        while (!vertex_set->empty()) {
            pool->submit_task([=]{
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
            });
            for (const auto &u: *vertex_set) {
                evicted_set->insert(u);
            }
            pool->barrier();
            auto location_vector = pool->split_task(vertex_set);
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (uint32_t i = 0; i < thread_number; ++i) {
               pool->submit_task([=] {
                    auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
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

    void branch_multiple_core_maintenance::remove_unsatisfied_vertices(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &vertex_set,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &candidate_map,
            const shared_ptr<unordered_set<uint32_t>> &evicted_set,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto global_mutex = make_shared<mutex>();
        auto thread_number = pool->get_thread_number();
        while (!vertex_set->empty()) {
            pool->submit_task([=]{
                for (const auto &u: *vertex_set) {
                    candidate_map->erase(u);
                }
            });
            for (const auto &u: *vertex_set) {
                evicted_set->insert(u);
            }
            pool->barrier();
            auto location_vector = pool->split_task(vertex_set);
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                        auto &u = *iter1;
                        auto &u_map = vertex_edge_index_map->at(u);
                        for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                            for (const auto &v: *iter2->second) {
                                if (candidate_map->count(v) && candidate_map->at(v) >=k) {

                                    vertex_mutex_map->at(v)->lock();
                                    --candidate_map->at(v);
                                    vertex_mutex_map->at(v)->unlock();

                                    if (candidate_map->at(v) < k) {
                                        sub_next_vertex_set->insert(v);
                                    }
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

    void branch_multiple_core_maintenance::init(const shared_ptr<temporal_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                const shared_ptr<thread_pool> &pool) {
        auto vertex_map = G->get_vertex_map();
        for (const auto &[u, u_vertex]: *vertex_map) {
            vertex_mutex_map->insert({u, shared_ptr<mutex>()});
            vertex_degree_map->insert({u, 0});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_map);

        for (uint32_t i = 0; i < thread_number; ++i) {
            pool->submit_task([=] {
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                    auto u = iter->first;
                    vertex_mutex_map->at(u) = make_shared<mutex>();
                }
            });
        }
        pool->barrier();
    }


    void branch_multiple_core_maintenance::init(const shared_ptr<temporal_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_size_map,
                                                const shared_ptr<thread_pool> &pool) {
        auto vertex_map = G->get_vertex_map();
        for (const auto &[u, u_vertex]: *vertex_map) {
            vertex_edge_size_map->insert({u, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_map);
        for (uint32_t i = 0; i < thread_number; ++i) {
            auto sub_begin = *location_vector->at(i);
            auto sub_end = *location_vector->at(i + 1);
            pool->submit_task([=] {
                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                    auto &[u, u_vertex] = *iter;
                    vertex_edge_size_map->at(u) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                    auto u_map = vertex_edge_size_map->at(u);
                    for (const auto &[v, e_set]: *u_vertex->get_neighbor_map()) {
                        auto h = e_set->size();
                        if (!u_map->count(h)) {
                            u_map->insert({h, make_shared<unordered_set<uint32_t>>()});
                        }
                        u_map->at(h)->insert(v);
                    }
                }
            });
        }
        pool->barrier();
    }

    void branch_multiple_core_maintenance::init(const shared_ptr<temporal_graph> &G,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_size_map,
                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                const shared_ptr<thread_pool> &pool) {
        auto vertex_map = G->get_vertex_map();
        for (const auto &[u, u_vertex]: *vertex_map) {
            vertex_mutex_map->insert({u, make_shared<mutex>()});
            vertex_edge_size_map->insert({u, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
            vertex_degree_map->insert({u, 0});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_map);
        for (uint32_t i = 0; i < thread_number; ++i) {
            auto sub_begin = *location_vector->at(i);
            auto sub_end = *location_vector->at(i + 1);
            pool->submit_task([=] {
                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                    auto &[u, u_vertex] = *iter;
                    vertex_mutex_map->at(u) = make_shared<mutex>();
                    vertex_edge_size_map->at(u) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                    auto u_map = vertex_edge_size_map->at(u);
                    for (const auto &[v, e_set]: *u_vertex->get_neighbor_map()) {
                        auto h = e_set->size();
                        if (!u_map->count(h)) {
                            u_map->insert({h, make_shared<unordered_set<uint32_t>>()});
                        }
                        u_map->at(h)->insert(v);
                    }
                }
            });
        }
        pool->barrier();
    }


    void branch_multiple_core_maintenance::insert(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto affected_set = make_shared<unordered_set<uint32_t>>();
        G->insert_edge_collection(edge_set);
        {
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();

                affected_set->insert(u);
                affected_set->insert(v);
            }

            for (const auto &u: *affected_set) {
                if (!vertex_index_map->count(u)) {
                    vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    vertex_index_map->at(u)->insert(1, 1);
                }
            }
        }

        auto previous_middle_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto max_delta = find_max_delta(vertex_index_map);

        auto global_mutex = make_shared<mutex>();
        for (uint32_t delta = 1; delta <= max_delta + 1; ++delta) {
            auto copy_previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(*previous_middle_candidate_map);
            pool->submit_task([=] {
                if (delta == 1) {
                    /**
                     * @brief update (1, h)-cores (h >= 2)
                     */
                    auto h_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                    for(const auto &u:*affected_set){
                        uint32_t max_h = 0;
                        for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                            if(affected_set->count(v) && v_edge_set->size() >= max_h){
                                max_h = v_edge_set->size();
                            }
                        }
                        if(max_h > 0){
                            if(!h_map->count(max_h)){
                                h_map->insert({max_h, make_shared<unordered_set<uint32_t>>()});
                            }
                            h_map->at(max_h)->insert(u);
                        }
                    }
                    global_mutex->lock();
                    for(const auto &[h, h_set]:*h_map){
                        for(const auto &u:*h_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                            new_vertex_index_map->at(u)->insert(1, h);
                        }
                    }
                    global_mutex->unlock();
                    /**
                     * @brief update (k, 1)-cores (k >= 2)
                     */
                    {
                        auto previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto max_k = find_max_k(affected_set, vertex_index_map, delta);

                        for (uint32_t k = 2; k <= max_k + 1; ++k) {
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                    *previous_candidate_map);

                            auto flag = left_candidate_graph(G, affected_set, vertex_index_map, vertex_degree_map,
                                                             previous_candidate_map,
                                                             candidate_map, k, 1);

                            if (flag) {

                                global_mutex->lock();
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(k, 1);
                                }
                                global_mutex->unlock();

                                swap(*previous_candidate_map, *candidate_map);
                            } else {
                                swap(*previous_candidate_map, *candidate_map);
                                break;
                            }
                        }
                        if (!previous_candidate_map->empty()) {
                            left_decomposition(G, global_mutex, previous_candidate_map, new_vertex_index_map, max_k + 2,
                                               1);
                        }
                    }
                }
                else {
                    /**
                     * @brief update (k, delta)-cores
                     */
                    {
                        auto previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                *copy_previous_candidate_map);
                        auto max_k = find_max_k(affected_set, vertex_index_map, delta);

                        for (uint32_t k = delta + 1; k <= max_k + 1; ++k) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                    *previous_candidate_map);
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto flag = left_candidate_graph(G, affected_set, vertex_index_map, vertex_degree_map,
                                                             previous_candidate_map,
                                                             candidate_map, k, delta);

                            if (flag) {
                                global_mutex->lock();
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(k, delta);
                                }
                                global_mutex->unlock();

                                swap(*previous_candidate_map, *candidate_map);
                            } else {
                                break;
                            }
                        }
                        if (!previous_candidate_map->empty()) {
                            left_decomposition(G, global_mutex, previous_candidate_map,
                                               new_vertex_index_map, max_k + 2,
                                               delta);
                        }
                    }
                    /**
                     * @brief update (delta, h)-cores
                     */
                    {
                        auto &previous_candidate_map = copy_previous_candidate_map;

                        auto max_h = find_max_h(affected_set, vertex_index_map, delta);

                        for (uint32_t h = delta + 1; h <= max_h + 1; ++h) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                    *previous_candidate_map);
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                            auto flag = right_candidate_graph(G, affected_set, vertex_index_map,
                                                              vertex_degree_map,
                                                              previous_candidate_map, candidate_map, delta, h);

                            if (flag) {
                                global_mutex->lock();
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(delta, h);
                                }
                                global_mutex->unlock();

                                swap(*previous_candidate_map, *candidate_map);
                            } else {
                                swap(*previous_candidate_map, *candidate_map);
                                break;
                            }
                        }

                        if (!previous_candidate_map->empty()) {
                            right_decomposition(G, global_mutex, previous_candidate_map,
                                                new_vertex_index_map, delta,
                                                max_h + 2);
                        }
                    }
                }
            });

            {
                auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*previous_middle_candidate_map);
                auto middle_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto flag = middle_candidate_graph(G, affected_set, vertex_index_map, vertex_degree_map, previous_middle_candidate_map,
                                                   middle_candidate_map, delta + 1, delta + 1);
                if (flag) {
                    global_mutex->lock();
                    for (const auto &[u, u_degree]: *middle_candidate_map) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                        }
                        new_vertex_index_map->at(u)->insert(delta + 1, delta + 1);
                    }
                    global_mutex->unlock();

                    swap(*previous_middle_candidate_map, *middle_candidate_map);
                } else {
                    swap(*previous_middle_candidate_map, *middle_candidate_map);
                    break;
                }
            }
        }

        if (!previous_middle_candidate_map->empty()) {
            middle_decomposition(G, global_mutex, previous_middle_candidate_map,
                                 new_vertex_index_map, max_delta + 2,
                                 pool);
        }
        pool->barrier();

        for (const auto &[u, u_index]: *new_vertex_index_map) {
            vertex_index_map->at(u)->merge_insert(u_index);
            u_index->clear();
        }
    }

    void branch_multiple_core_maintenance::insert(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t,
                                                          shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto affected_set = make_shared<unordered_set<uint32_t>>();
        update_vertex_edge_index_for_insertion(G, edge_set, vertex_edge_index_map, affected_set);
        {
            for (const auto &u: *affected_set) {
                if (!vertex_index_map->count(u)) {
                    vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    vertex_index_map->at(u)->insert(1, 1);
                }
            }
        }

        auto previous_middle_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto max_delta = find_max_delta(vertex_index_map);

        auto global_mutex = make_shared<mutex>();
        for (uint32_t delta = 1; delta <= max_delta + 1; ++delta) {
            auto copy_previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(*previous_middle_candidate_map);
            pool->submit_task([=] {
                if (delta == 1) {
                    /**
                     * @brief update (1, h)-cores (h >= 2)
                     */
                    auto h_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                    for(const auto &u:*affected_set){
                        auto max_h = vertex_edge_index_map->at(u)->rbegin()->first;
                        if(max_h > 0){
                            if(!h_map->count(max_h)){
                                h_map->insert({max_h, make_shared<unordered_set<uint32_t>>()});
                            }
                            h_map->at(max_h)->insert(u);
                        }
                    }
                    global_mutex->lock();
                    for(const auto &[h, h_set]:*h_map){
                        for(const auto &u:*h_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                            new_vertex_index_map->at(u)->insert(1, h);
                        }
                    }
                    global_mutex->unlock();
                    /**
                     * @brief update (k, 1)-cores (k >= 2)
                     */
                    {
                        auto previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto max_k = find_max_k(affected_set, vertex_index_map, delta);

                        for (uint32_t k = 2; k <= max_k + 1; ++k) {
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                    *previous_candidate_map);

                            auto flag = left_candidate_graph(vertex_edge_index_map, affected_set, vertex_index_map, vertex_degree_map,
                                                             previous_candidate_map,
                                                             candidate_map, k, 1);

                            if (flag) {

                                global_mutex->lock();
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(k, 1);
                                }
                                global_mutex->unlock();

                                swap(*previous_candidate_map, *candidate_map);
                            } else {
                                swap(*previous_candidate_map, *candidate_map);
                                break;
                            }
                        }
                        if (!previous_candidate_map->empty()) {
                            left_decomposition(vertex_edge_index_map, global_mutex, previous_candidate_map, new_vertex_index_map, max_k + 2,
                                               1);
                        }
                    }
                }
                else {
                    /**
                     * @brief update (k, delta)-cores
                     */
                    {
                        auto previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                *copy_previous_candidate_map);
                        auto max_k = find_max_k(affected_set, vertex_index_map, delta);

                        for (uint32_t k = delta + 1; k <= max_k + 1; ++k) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                    *previous_candidate_map);
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto flag = left_candidate_graph(vertex_edge_index_map, affected_set, vertex_index_map, vertex_degree_map,
                                                             previous_candidate_map,
                                                             candidate_map, k, delta);

                            if (flag) {
                                global_mutex->lock();
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(k, delta);
                                }
                                global_mutex->unlock();

                                swap(*previous_candidate_map, *candidate_map);
                            } else {
                                swap(*previous_candidate_map, *candidate_map);
                                break;
                            }
                        }
                        if (!previous_candidate_map->empty()) {
                            left_decomposition(vertex_edge_index_map, global_mutex, previous_candidate_map,
                                               new_vertex_index_map, max_k + 2,
                                               delta);
                        }
                    }
                    /**
                     * @brief update (delta, h)-cores
                     */
                    {
                        auto &previous_candidate_map = copy_previous_candidate_map;

                        auto max_h = find_max_h(affected_set, vertex_index_map, delta);

                        for (uint32_t h = delta + 1; h <= max_h + 1; ++h) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                    *previous_candidate_map);
                            auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                            auto flag = right_candidate_graph(vertex_edge_index_map, affected_set, vertex_index_map,
                                                              vertex_degree_map,
                                                              previous_candidate_map, candidate_map, delta, h);

                            if (flag) {
                                global_mutex->lock();
                                for (const auto &[u, u_degree]: *candidate_map) {
                                    if(!new_vertex_index_map->count(u)){
                                        new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                    }
                                    new_vertex_index_map->at(u)->insert(delta, h);
                                }
                                global_mutex->unlock();

                                swap(*previous_candidate_map, *candidate_map);
                            } else {
                                swap(*previous_candidate_map, *candidate_map);
                                break;
                            }
                        }

                        if (!previous_candidate_map->empty()) {
                            right_decomposition(vertex_edge_index_map, global_mutex, previous_candidate_map,
                                                new_vertex_index_map, delta,
                                                max_h + 2);
                        }
                    }
                }
            });

            {
                auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*previous_middle_candidate_map);
                auto middle_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto flag = middle_candidate_graph(vertex_edge_index_map, affected_set, vertex_index_map, vertex_degree_map, previous_middle_candidate_map,
                                                   middle_candidate_map, delta + 1, delta + 1);
                if (flag) {
                    global_mutex->lock();
                    for (const auto &[u, u_degree]: *middle_candidate_map) {
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                        }
                        new_vertex_index_map->at(u)->insert(delta + 1, delta + 1);
                    }
                    global_mutex->unlock();

                    swap(*previous_middle_candidate_map, *middle_candidate_map);
                } else {
                    swap(*previous_middle_candidate_map, *middle_candidate_map);
                    break;
                }
            }
        }

        if (!previous_middle_candidate_map->empty()) {
            middle_decomposition(vertex_edge_index_map, global_mutex, previous_middle_candidate_map,
                                 new_vertex_index_map, max_delta + 2,
                                 pool);
        }
        pool->barrier();

        for (const auto &[u, u_index]: *new_vertex_index_map) {
            vertex_index_map->at(u)->merge_insert(u_index);
            u_index->clear();
        }
    }

    void branch_multiple_core_maintenance::insert(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> & vertex_degree_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto affected_set = make_shared<unordered_set<uint32_t>>();
        pool->submit_task([=]{
            G->insert_edge_collection(edge_set);
        });

        {
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();

                affected_set->insert(u);
                affected_set->insert(v);
            }

            for (const auto &u: *affected_set) {
                if (!vertex_index_map->count(u)) {
                    vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    vertex_index_map->at(u)->insert(1, 1);

                    vertex_mutex_map->insert({u, make_shared<mutex>()});
                    vertex_degree_map->insert({u, 0});
                }
            }
        }
        pool->barrier();

        auto previous_middle_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto max_delta = find_max_delta(vertex_index_map, pool);

        for (uint32_t delta = 1; delta <= max_delta + 1; ++delta) {
            if (delta == 1) {
                /**
                 * @brief update (1, h)-cores (h >= 2)
                 */
                for(const auto &u:*affected_set){
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                auto thread_number = pool->get_thread_number();
                auto location_vector = pool->split_task(affected_set);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &u = *iter;
                            if(!new_vertex_index_map->at(u)){
                                new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                            }

                            uint32_t max_h = 0;
                            for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                                if(affected_set->count(v) && v_edge_set->size() >= max_h){
                                    max_h = v_edge_set->size();
                                }
                            }
                            if (max_h > 0) {
                                new_vertex_index_map->at(u)->insert(1, max_h);
                            }
                        }
                    });
                }
                pool->barrier();
                /**
                 * @brief update (k, 1)-cores (k >= 2)
                 */
                {
                    auto previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto max_k = find_max_k(affected_set, vertex_index_map, delta, pool);

                    for (uint32_t k = 2; k <= max_k + 1; ++k) {
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto flag = left_candidate_graph(G, vertex_mutex_map, affected_set,
                                                         vertex_index_map,
                                                         vertex_degree_map,
                                                         previous_candidate_map,
                                                         candidate_map, k, 1, pool);

                        if (flag) {
                            for(const auto &[u, u_degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                                }
                            }
                            insertion_assign(candidate_map, new_vertex_index_map, k, 1, pool);
                            swap(*previous_candidate_map, *candidate_map);
                        } else {
                            swap(*previous_candidate_map, *candidate_map);
                            break;
                        }
                    }
                    if (!previous_candidate_map->empty()) {
                        left_decomposition(G, vertex_mutex_map, previous_candidate_map,
                                           new_vertex_index_map, max_k + 2,
                                           1, pool);
                    }
                }
            } else {
                auto previous_candidate_map1 = make_shared<unordered_map<uint32_t, uint32_t>>();
                pool->submit_task([=]{
                    copy(previous_middle_candidate_map->begin(), previous_middle_candidate_map->end(), std::inserter(*previous_candidate_map1, previous_candidate_map1->end()));
                });
                auto previous_candidate_map2 = make_shared<unordered_map<uint32_t, uint32_t>>(*previous_middle_candidate_map);
                pool->barrier();
                /**
                 * @brief update (k, delta)-cores
                 */
                {

                    auto &previous_candidate_map = previous_candidate_map1;
                    auto max_k = find_max_k(affected_set, vertex_index_map, delta, pool);
                    for (uint32_t k = delta + 1; k <= max_k + 1; ++k) {
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto flag = left_candidate_graph(G, vertex_mutex_map, affected_set,
                                                         vertex_index_map,
                                                         vertex_degree_map,
                                                         previous_candidate_map,
                                                         candidate_map, k, delta, pool);

                        if (flag) {
                            for(const auto &[u, u_degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                                }
                            }
                            insertion_assign(candidate_map, new_vertex_index_map, k, delta, pool);
                            swap(*previous_candidate_map, *candidate_map);
                        } else {
                            swap(*previous_candidate_map, *candidate_map);
                            break;
                        }
                    }
                    if (!previous_candidate_map->empty()) {
                        left_decomposition(G, vertex_mutex_map, previous_candidate_map,
                                           new_vertex_index_map, max_k + 2,
                                           delta, pool);
                    }
                }
                /**
                 * @brief update (delta, h)-cores
                 */
                {
                    auto &previous_candidate_map = previous_candidate_map2;
                    auto max_h = find_max_h(affected_set, vertex_index_map, delta, pool);
                    for (uint32_t h = delta + 1; h <= max_h + 1; ++h) {
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                        auto flag = right_candidate_graph(G, vertex_mutex_map, affected_set,
                                                          vertex_index_map,
                                                          vertex_degree_map,
                                                          previous_candidate_map, candidate_map, delta, h, pool);

                        if (flag) {
                            for(const auto &[u, u_degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                                }
                            }
                            insertion_assign(candidate_map, new_vertex_index_map, delta, h, pool);
                            swap(*previous_candidate_map, *candidate_map);
                        } else {
                            swap(*previous_candidate_map, *candidate_map);
                            break;
                        }
                    }

                    if (!previous_candidate_map->empty()) {
                        right_decomposition(G, vertex_mutex_map, previous_candidate_map,
                                            new_vertex_index_map, delta,
                                            max_h + 2, pool);
                    }
                }
            }

            {
                auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto flag = middle_candidate_graph(G, vertex_mutex_map, affected_set,
                                                   vertex_index_map,
                                                   vertex_degree_map, previous_middle_candidate_map,
                                                   candidate_map, delta + 1, delta + 1, pool);
                if (flag) {
                    for(const auto &[u, u_degree]:*candidate_map){
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }

                    insertion_assign(candidate_map, new_vertex_index_map, delta + 1, delta + 1, pool);
                    swap(*previous_middle_candidate_map, *candidate_map);
                } else {
                    swap(*previous_middle_candidate_map, *candidate_map);
                    break;
                }
            }
        }

        if (!previous_middle_candidate_map->empty()) {
            middle_decomposition(G, vertex_mutex_map, previous_middle_candidate_map,
                                 new_vertex_index_map, max_delta + 2,
                                 pool);
        }

        insertion_merge(new_vertex_index_map, vertex_index_map, pool);
    }

    void branch_multiple_core_maintenance::insert(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t,
                                                          shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> & vertex_degree_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto affected_set = make_shared<unordered_set<uint32_t>>();
        update_vertex_edge_index_for_insertion(G, edge_set, vertex_edge_index_map, affected_set, pool);
        {
            for (const auto &u: *affected_set) {
                if (!vertex_index_map->count(u)) {
                    vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                    vertex_index_map->at(u)->insert(1, 1);

                    vertex_mutex_map->insert({u, make_shared<mutex>()});
                    vertex_degree_map->insert({u, 0});
                }
            }
        }

        auto previous_middle_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto max_delta = find_max_delta(vertex_index_map, pool);

        for (uint32_t delta = 1; delta <= max_delta + 1; ++delta) {
            if (delta == 1) {
                /**
                 * @brief update (1, h)-cores (h >= 2)
                 */
                for(const auto &u:*affected_set){
                    if(!new_vertex_index_map->count(u)){
                        new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                    }
                }
                auto thread_number = pool->get_thread_number();
                auto location_vector = pool->split_task(affected_set);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &u = *iter;
                            if(!new_vertex_index_map->at(u)){
                                new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                            }

                            auto max_h = vertex_edge_index_map->at(u)->rbegin()->first;
                            if (max_h > 0) {
                                new_vertex_index_map->at(u)->insert(1, max_h);
                            }
                        }
                    });
                }
                pool->barrier();
                /**
                 * @brief update (k, 1)-cores (k >= 2)
                 */
                {
                    auto previous_candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto max_k = find_max_k(affected_set, vertex_index_map, delta, pool);

                    for (uint32_t k = 2; k <= max_k + 1; ++k) {
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto flag = left_candidate_graph(vertex_edge_index_map, vertex_mutex_map, affected_set, vertex_index_map,
                                                         vertex_degree_map,
                                                         previous_candidate_map,
                                                         candidate_map, k, 1, pool);

                        if (flag) {
                            for(const auto &[u, u_degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                                }
                            }

                            insertion_assign(candidate_map, new_vertex_index_map, k, 1, pool);
                            swap(*previous_candidate_map, *candidate_map);
                        } else {
                            swap(*previous_candidate_map, *candidate_map);
                            break;
                        }
                    }
                    if (!previous_candidate_map->empty()) {
                        left_decomposition(vertex_edge_index_map, vertex_mutex_map, previous_candidate_map, new_vertex_index_map, max_k + 2,
                                           1, pool);
                    }
                }
            } else {
                /**
                 * @brief update (k, delta)-cores
                 */
                auto previous_candidate_map1 = make_shared<unordered_map<uint32_t, uint32_t>>(*previous_middle_candidate_map);
                pool->submit_task([=]{
                    copy(previous_middle_candidate_map->begin(), previous_middle_candidate_map->end(), inserter(*previous_candidate_map1, previous_candidate_map1->end()));
                });
                auto previous_candidate_map2 = make_shared<unordered_map<uint32_t, uint32_t>>(*previous_middle_candidate_map);
                pool->barrier();
                {
                    auto &previous_candidate_map = previous_candidate_map1;

                    auto max_k = find_max_k(affected_set, vertex_index_map, delta, pool);
                    for (uint32_t k = delta + 1; k <= max_k + 1; ++k) {
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto flag = left_candidate_graph(vertex_edge_index_map, vertex_mutex_map, affected_set, vertex_index_map,
                                                         vertex_degree_map,
                                                         previous_candidate_map,
                                                         candidate_map, k, delta, pool);

                        if (flag) {
                            for(const auto &[u, u_degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                                }
                            }

                            insertion_assign(candidate_map, new_vertex_index_map, k, delta, pool);
                            swap(*previous_candidate_map, *candidate_map);
                        } else {
                            swap(*previous_candidate_map, *candidate_map);
                            break;
                        }
                    }
                    if (!previous_candidate_map->empty()) {
                        left_decomposition(vertex_edge_index_map, vertex_mutex_map, previous_candidate_map,
                                           new_vertex_index_map, max_k + 2,
                                           delta, pool);
                    }
                }
                /**
                 * @brief update (delta, h)-cores
                 */
                {
                    auto &previous_candidate_map = previous_candidate_map2;

                    auto max_h = find_max_h(affected_set, vertex_index_map, delta, pool);
                    for (uint32_t h = delta + 1; h <= max_h + 1; ++h) {
                        auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                        auto flag = right_candidate_graph(vertex_edge_index_map, vertex_mutex_map, affected_set, vertex_index_map,
                                                          vertex_degree_map,
                                                          previous_candidate_map, candidate_map, delta, h, pool);

                        if (flag) {
                            for(const auto &[u, u_degree]:*candidate_map){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                                }
                            }

                            insertion_assign(candidate_map, new_vertex_index_map, delta,  h, pool);
                            swap(*previous_candidate_map, *candidate_map);
                        } else {
                            swap(*previous_candidate_map, *candidate_map);
                            break;
                        }
                    }

                    if (!previous_candidate_map->empty()) {
                        right_decomposition(vertex_edge_index_map, vertex_mutex_map, previous_candidate_map,
                                            new_vertex_index_map, delta,
                                            max_h + 2, pool);
                    }
                }
            }

            {
                auto candidate_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto flag = middle_candidate_graph(vertex_edge_index_map, vertex_mutex_map, affected_set, vertex_index_map,
                                                   vertex_degree_map, previous_middle_candidate_map,
                                                   candidate_map, delta + 1, delta + 1, pool);
                if (flag) {
                    for(const auto &[u, u_degree]:*candidate_map){
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }

                    insertion_assign(candidate_map, new_vertex_index_map, delta + 1, delta + 1, pool);
                    swap(*previous_middle_candidate_map, *candidate_map);
                } else {
                    swap(*previous_middle_candidate_map, *candidate_map);
                    break;
                }
            }
        }

        if (!previous_middle_candidate_map->empty()) {
            middle_decomposition(vertex_edge_index_map, vertex_mutex_map, previous_middle_candidate_map,
                                 new_vertex_index_map, max_delta + 2,
                                 pool);
        }

        insertion_merge(new_vertex_index_map, vertex_index_map, pool);
    }

    void branch_multiple_core_maintenance::remove(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();
        auto affected_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        {
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();

                affected_set->insert(u);
                affected_set->insert(v);
            }

            G->remove_edge_collection(edge_set, isolated_vertex_set);

            for (const auto &u: *isolated_vertex_set) {
                affected_set->erase(u);
            }
        }

        auto max_delta = find_max_delta(affected_set, vertex_index_map);

        auto global_mutex = make_shared<mutex>();
        for (uint32_t delta = max_delta; delta >= 1; --delta) {
            pool->submit_task([=] {
                if (delta == 1) {
                    {
                        /**
                         * @brief update (1, h)-cores
                         */
                        auto h_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                        for (const auto &u: *affected_set) {
                            uint32_t max_h = 0;
                            for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                auto h = v_edge_set->size();
                                if (h > max_h) {
                                    max_h = v_edge_set->size();
                                }
                            }
                            if(max_h > 0){
                                if(!h_map->count(max_h)){
                                    h_map->insert({max_h, make_shared<unordered_set<uint32_t>>()});
                                }
                                h_map->at(max_h)->insert(u);
                            }
                        }
                        global_mutex->lock();
                        for(const auto&[h, h_set]:*h_map){
                            for(const auto&u:*h_set){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(1, h + 1);
                            }
                        }
                        global_mutex->unlock();
                    }
                    /**
                     * @brief update (k, 1)-cores (k >=2)
                     */
                    {
                        auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                        auto max_k = find_max_k(affected_set, vertex_index_map, 1);
                        for (uint32_t  k = max_k; k > 1; --k) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto removed_set = left_removal_core(G, affected_set, vertex_index_map, vertex_degree_map,
                                                                 previous_removed_set, k, 1);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(k, 1);
                            }
                            global_mutex->unlock();

                            swap(*previous_removed_set, *removed_set);
                        }
                    }
                } else {
                    {
                        auto max_k = find_max_k(affected_set, vertex_index_map, delta);

                        auto previous_removed_set = make_shared<unordered_set<uint32_t>>();

                        for (uint32_t k = max_k; k >= delta; --k) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto removed_set = left_removal_core(G, affected_set, vertex_index_map,
                                                                 vertex_degree_map,
                                                                 previous_removed_set, k,
                                                                 delta);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(k, delta);
                            }
                            global_mutex->unlock();

                            swap(*removed_set, *previous_removed_set);
                        }
                    }

                    {
                        auto max_h = find_max_h(affected_set, vertex_index_map, delta);

                        auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                        for (uint32_t h = max_h; h > delta; --h) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto removed_set = right_removal_core(G, affected_set, vertex_index_map,
                                                                  vertex_degree_map,
                                                                  previous_removed_set, delta, h);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(delta, h);
                            }
                            global_mutex->unlock();

                            swap(*removed_set, *previous_removed_set);
                        }
                    }
                }
            });
        }
        pool->barrier();

        for (const auto &[u, u_index]: *new_vertex_index_map) {
            vertex_index_map->at(u)->merge_remove(u_index);
        }

        for (const auto &u: *isolated_vertex_set) {
            vertex_index_map->erase(u);
        }
    }

    void branch_multiple_core_maintenance::remove(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto affected_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        update_vertex_edge_index_for_removal(G, edge_set, vertex_edge_index_map, affected_set, isolated_vertex_set);

        auto max_delta = find_max_delta(affected_set, vertex_index_map);

        auto global_mutex = make_shared<mutex>();
        auto previous_middle_removed_set = make_shared<unordered_set<uint32_t>>();
        for (uint32_t delta = max_delta; delta >= 1; --delta) {
            pool->submit_task([=] {
                if (delta == 1) {
                    {
                        /**
                         * @brief update (1, h)-cores
                         */
                        auto h_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                        for (const auto &u: *affected_set) {
                            auto max_h = vertex_edge_index_map->at(u)->rbegin()->first;
                            if(max_h > 0){
                                if(!h_map->count(max_h)){
                                    h_map->insert({max_h, make_shared<unordered_set<uint32_t>>()});
                                }
                                h_map->at(max_h)->insert(u);
                            }
                        }
                        global_mutex->lock();
                        for(const auto&[h, h_set]:*h_map){
                            for(const auto&u:*h_set){
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(1, h + 1);
                            }
                        }
                        global_mutex->unlock();
                    }
                    /**
                     * @brief update (k, 1)-cores (k >=2)
                     */
                    {
                        auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                        auto max_k = find_max_k(affected_set, vertex_index_map, 1);
                        for (uint32_t  k = max_k; k > 1; --k) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto removed_set = left_removal_core(vertex_edge_index_map, affected_set, vertex_index_map, vertex_degree_map,
                                                                 previous_removed_set, k, 1);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(k, 1);
                            }
                            global_mutex->unlock();

                            swap(*previous_removed_set, *removed_set);
                        }
                    }
                } else {
                    {
                        auto max_k = find_max_k(affected_set, vertex_index_map, delta);

                        auto previous_removed_set = make_shared<unordered_set<uint32_t>>();

                        for (uint32_t k = max_k; k >= delta; --k) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto removed_set = left_removal_core(vertex_edge_index_map, affected_set, vertex_index_map,
                                                                 vertex_degree_map,
                                                                 previous_removed_set,
                                                                 k,
                                                                 delta);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(k, delta);
                            }
                            global_mutex->unlock();

                            swap(*removed_set, *previous_removed_set);
                        }
                    }

                    {
                        auto max_h = find_max_h(affected_set, vertex_index_map, delta);

                        auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                        for (uint32_t h = max_h; h > delta; --h) {
                            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                            auto removed_set = right_removal_core(vertex_edge_index_map, affected_set, vertex_index_map,
                                                                  vertex_degree_map,
                                                                  previous_removed_set, delta, h);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                if(!new_vertex_index_map->count(u)){
                                    new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                                }
                                new_vertex_index_map->at(u)->remove(delta, h);
                            }
                            global_mutex->unlock();

                            swap(*removed_set, *previous_removed_set);
                        }
                    }
                }
            });
        }
        pool->barrier();

        for (const auto &[u, u_index]: *new_vertex_index_map) {
            vertex_index_map->at(u)->merge_remove(u_index);
        }

        for (const auto &u: *isolated_vertex_set) {
            vertex_index_map->erase(u);
            vertex_edge_index_map->erase(u);
        }
    }

    void branch_multiple_core_maintenance::remove(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> & vertex_degree_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto affected_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        {
            pool->submit_task([=]{
                G->remove_edge_collection(edge_set, isolated_vertex_set);
            });
            for (const auto &e: *edge_set) {
                auto u = e->get_source_vertex_id();
                auto v = e->get_destination_vertex_id();

                affected_set->insert(u);
                affected_set->insert(v);
            }
            pool->barrier();


            for (const auto &u: *isolated_vertex_set) {
                affected_set->erase(u);
            }
        }

        auto max_delta = find_max_delta(affected_set, vertex_index_map, pool);

        auto previous_middle_removed_set = make_shared<unordered_set<uint32_t>>();
        for (uint32_t delta = max_delta; delta >= 1; --delta) {
            if (delta == 1) {
                {
                    /**
                     * @brief update (1, h)-cores
                     */
                    for(const auto &u:*affected_set){
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                        }
                    }

                    auto thread_number = pool->get_thread_number();
                    auto location_vector = pool->split_task(affected_set);
                    for (uint32_t i = 0; i < thread_number; ++i) {
                        pool->submit_task([=] {
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                auto &u = *iter;
                                if(!new_vertex_index_map->at(u)){
                                    new_vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                                }
                                uint32_t max_h = 0;
                                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                    auto h = v_edge_set->size();
                                    if (h > max_h) {
                                        max_h = h;
                                    }
                                }
                                new_vertex_index_map->at(u)->remove(1, max_h + 1);
                            }
                        });
                    }
                    pool->barrier();
                }
                /**
                 * @brief update (k, 1)-cores (k >=2)
                 */
                {
                    auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                    auto max_k = find_max_k(affected_set, vertex_index_map, 1, pool);
                    for (uint32_t k = max_k; k > 1; --k) {
                        auto removed_set = left_removal_core(G, vertex_mutex_map, affected_set, vertex_index_map,
                                                             vertex_degree_map,
                                                             previous_removed_set, k, 1, pool);

                        for(const auto &u:*removed_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                            }
                        }

                        removal_assign(removed_set, new_vertex_index_map, k, 1, pool);
                        swap(*previous_removed_set, *removed_set);
                    }
                }
            } else {
                {
                    auto max_k = find_max_k(affected_set, vertex_index_map, delta, pool);
                    auto previous_removed_set = make_shared<unordered_set<uint32_t>>();

                    for (uint32_t k = max_k; k >= delta; --k) {
                        auto removed_set = left_removal_core(G, vertex_mutex_map, affected_set, vertex_index_map,
                                                             vertex_degree_map,
                                                             previous_removed_set, k,
                                                             delta, pool);

                        for(const auto &u:*removed_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                            }
                        }

                        removal_assign(removed_set, new_vertex_index_map, k, delta, pool);
                        swap(*removed_set, *previous_removed_set);
                    }
                }

                {
                    auto max_h = find_max_h(affected_set, vertex_index_map, delta, pool);

                    auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                    for (uint32_t h = max_h; h > delta; --h) {
                        auto removed_set = right_removal_core(G, vertex_mutex_map, affected_set, vertex_index_map,
                                                              vertex_degree_map,
                                                              previous_removed_set, delta, h, pool);

                        for(const auto &u:*removed_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
                            }
                        }

                        removal_assign(removed_set, new_vertex_index_map, delta, h, pool);
                        swap(*removed_set, *previous_removed_set);
                    }
                }
            }
        }

        removal_merge(new_vertex_index_map, vertex_index_map, pool);

        for (const auto &u: *isolated_vertex_set) {
            vertex_index_map->erase(u);
            vertex_mutex_map->erase(u);
            vertex_degree_map->erase(u);
        }
    }

    void branch_multiple_core_maintenance::remove(const shared_ptr<temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                  const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                  const shared_ptr<thread_pool> &pool) {
        auto new_vertex_index_map =  make_shared<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>();

        auto affected_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        update_vertex_edge_index_for_removal(G, edge_set, vertex_edge_index_map, affected_set, isolated_vertex_set, pool);

        auto max_delta = find_max_delta(affected_set, vertex_index_map, pool);

        auto previous_middle_removed_set = make_shared<unordered_set<uint32_t>>();
        for (uint32_t delta = max_delta; delta >= 1; --delta) {
            if (delta == 1) {
                {
                    /**
                     * @brief update (1, h)-cores
                     */
                    for(const auto &u:*affected_set){
                        if(!new_vertex_index_map->count(u)){
                            new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                        }
                    }

                    auto thread_number = pool->get_thread_number();
                    auto location_vector = pool->split_task(affected_set);
                    for(uint32_t i = 0; i < thread_number; ++i){
                        pool->submit_task([=]{
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &u = *iter;
                                auto max_h = vertex_edge_index_map->at(u)->rbegin()->first;
                                new_vertex_index_map->at(u)->remove(1, max_h + 1);
                            }
                        });
                    }
                    pool->barrier();
                }
                /**
                 * @brief update (k, 1)-cores (k >=2)
                 */
                {
                    auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                    auto max_k = find_max_k(affected_set, vertex_index_map, 1, pool);
                    for (uint32_t  k = max_k; k > 1; --k) {
                        auto removed_set = left_removal_core(vertex_edge_index_map, vertex_mutex_map, affected_set, vertex_index_map, vertex_degree_map,
                                                             previous_removed_set, k, 1, pool);

                        for(const auto &u:*removed_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                        }

                        removal_assign(removed_set, new_vertex_index_map, k, 1, pool);
                        swap(*previous_removed_set, *removed_set);
                    }
                }
            } else {
                {
                    auto max_k = find_max_k(affected_set, vertex_index_map, delta, pool);

                    auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                    for (uint32_t k = max_k; k >= delta; --k) {
                        auto removed_set = left_removal_core(vertex_edge_index_map, vertex_mutex_map, affected_set, vertex_index_map,
                                                             vertex_degree_map,
                                                             previous_removed_set, k,
                                                             delta, pool);

                        for(const auto &u:*removed_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                        }

                        removal_assign(removed_set, new_vertex_index_map, k, delta, pool);
                        swap(*removed_set, *previous_removed_set);
                    }
                }

                {
                    auto max_h = find_max_h(affected_set, vertex_index_map, delta, pool);

                    auto previous_removed_set = make_shared<unordered_set<uint32_t>>();
                    for (uint32_t h = max_h; h > delta; --h) {
                        auto removed_set = right_removal_core(vertex_edge_index_map, vertex_mutex_map, affected_set, vertex_index_map,
                                                              vertex_degree_map,
                                                              previous_removed_set, delta, h, pool);

                        for(const auto &u:*removed_set){
                            if(!new_vertex_index_map->count(u)){
                                new_vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
                            }
                        }

                        removal_assign(removed_set, new_vertex_index_map, delta, h, pool);
                        swap(*removed_set, *previous_removed_set);
                    }
                }
            }
        }

        removal_merge(new_vertex_index_map, vertex_index_map, pool);

        for (const auto &u: *isolated_vertex_set) {
            vertex_index_map->erase(u);
            vertex_edge_index_map->erase(u);
            vertex_mutex_map->erase(u);
            vertex_degree_map->erase(u);
        }
    }

    shared_ptr<unordered_set<uint32_t>>
    branch_multiple_core_maintenance::left_removal_core(const shared_ptr<temporal_graph> &G,
                                                        const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                        const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                        uint32_t k,
                                                        uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &u: *affected_set) {
            if (vertex_index_map->at(u)->equal(k, h) || previous_removed_set->count(u)) {
                vertex_degree_map->insert({u, compute_core_degree(G, vertex_index_map, u, k, h)});

                if (vertex_degree_map->at(u) < k) {
                    vertex_set->insert(u);
                }
            }
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                    for (const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()) {
                        if ((vertex_index_map->at(v)->equal(k, h) || previous_removed_set->count(v)) && v_edge_set->size() >= h
                            && !vertex_set->count(v) && !removed_set->count(v)) {
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
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_maintenance::left_removal_core(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &u: *affected_set) {
            if (vertex_index_map->at(u)->equal(k, h) || previous_removed_set->count(u)) {
                vertex_degree_map->insert({u, compute_core_degree(vertex_edge_index_map, vertex_index_map, u, k, h)});

                if (vertex_degree_map->at(u) < k) {
                    vertex_set->insert(u);
                }
            }
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                auto edge_map = vertex_edge_index_map->at(u);
                for (auto iter = edge_map->lower_bound(h); iter != edge_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if ((vertex_index_map->at(v)->equal(k, h) || previous_removed_set->count(v))
                            && !vertex_set->count(v) && !removed_set->count(v)) {
                            if (!vertex_degree_map->count(v)) {
                                vertex_degree_map->insert({v, compute_core_degree(vertex_edge_index_map, vertex_index_map, v, k, h)});
                            }

                            --vertex_degree_map->at(v);
                            if (vertex_degree_map->at(v) < k) {
                                next_vertex_set->insert(v);
                            }
                        }
                    }
                }
            }
            removed_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>>
    branch_multiple_core_maintenance::left_removal_core(const shared_ptr<temporal_graph> &G,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                        const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                        const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                        uint32_t k,
                                                        uint32_t h,
                                                        const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->equal(k, h) || previous_removed_set->count(u)) {
                            sub_clear_set->insert(u);
                            vertex_degree_map->at(u) = compute_core_degree(G, vertex_index_map, u, k, h);

                            if (vertex_degree_map->at(u) < k) {
                                sub_vertex_set->insert(u);
                            }
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        {
            while (!vertex_set->empty()) {
                auto location_vector = pool->split_task(vertex_set);
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                        auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);
                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &u = *iter;
                            for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                if ((vertex_index_map->at(v)->equal(k, h) || previous_removed_set->count(v))
                                    && e_set->size() >= h && !vertex_set->count(v) && !removed_set->count(v)) {
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
        }
        assign(clear_set, vertex_degree_map, 0, pool);
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_maintenance::left_removal_core(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->equal(k, h) || previous_removed_set->count(u)) {
                            sub_clear_set->insert(u);
                            vertex_degree_map->at(u) = compute_core_degree(vertex_edge_index_map, vertex_index_map, u, k, h);

                            if (vertex_degree_map->at(u) < k) {
                                sub_vertex_set->insert(u);
                            }
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        {
            while (!vertex_set->empty()) {
                auto location_vector = pool->split_task(vertex_set);
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                        auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);
                        for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                            auto &u = *iter1;
                            auto &u_map = vertex_edge_index_map->at(u);
                            for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                                for (const auto &v: *iter2->second) {
                                    if ((vertex_index_map->at(v)->equal(k, h) || previous_removed_set->count(v))
                                        && !vertex_set->count(v) && !removed_set->count(v)) {

                                        vertex_mutex_map->at(v)->lock();
                                        if (vertex_degree_map->at(v) == 0) {
                                            sub_clear_set->insert(v);
                                            vertex_degree_map->at(v) = compute_core_degree(vertex_edge_index_map,
                                                                                           vertex_index_map, v, k, h);
                                        }
                                        --vertex_degree_map->at(v);
                                        vertex_mutex_map->at(v)->unlock();

                                        if (vertex_degree_map->at(v) < k) {
                                            sub_next_vertex_set->insert(v);
                                        }
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
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return removed_set;
    }


    shared_ptr<unordered_set<uint32_t>>
    branch_multiple_core_maintenance::middle_removal_core(const shared_ptr<temporal_graph> &G,
                                                          const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                          const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                          uint32_t k,
                                                          uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &u: *affected_set) {
            if (vertex_index_map->at(u)->count(k, h) && (!vertex_index_map->at(u)->count(k + 1, h + 1) || previous_removed_set->count(u))) {
                vertex_degree_map->insert({u, compute_core_degree(G, vertex_index_map, u, k, h)});

                if (vertex_degree_map->at(u) < k) {
                    vertex_set->insert(u);
                }
            }
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_index_map->at(v)->count(k, h) &&
                    (!vertex_index_map->at(v)->count(k + 1, h + 1) || previous_removed_set->count(v))
                        && e_set->size() >= h && !vertex_set->count(v) && !removed_set->count(v)) {
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
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_maintenance::middle_removal_core(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &u: *affected_set) {
            if (vertex_index_map->at(u)->count(k, h) && (!vertex_index_map->at(u)->count(k + 1, h + 1) || previous_removed_set->count(u))) {
                vertex_degree_map->insert({u, compute_core_degree(vertex_edge_index_map, vertex_index_map, u, k, h)});

                if (vertex_degree_map->at(u) < k) {
                    vertex_set->insert(u);
                }
            }
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                auto edge_map = vertex_edge_index_map->at(u);
                for (auto iter = edge_map->lower_bound(h); iter != edge_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if (vertex_index_map->at(v)->count(k, h) &&(!vertex_index_map->at(v)->count(k + 1, h + 1) || previous_removed_set->count(v))
                            && !vertex_set->count(v) && !removed_set->count(v)) {
                            if (!vertex_degree_map->count(v)) {
                                vertex_degree_map->insert({v, compute_core_degree(vertex_edge_index_map, vertex_index_map, v, k, h)});
                            }

                            --vertex_degree_map->at(v);
                            if (vertex_degree_map->at(v) < k) {
                                next_vertex_set->insert(v);
                            }
                        }
                    }
                }
            }
            removed_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>>
    branch_multiple_core_maintenance::middle_removal_core(const shared_ptr<temporal_graph> &G,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                          const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                          const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                          uint32_t k,
                                                          uint32_t h,
                                                          const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->count(k, h) && (!vertex_index_map->at(u)->count(k + 1, h + 1) || previous_removed_set->count(u))) {
                            sub_clear_set->insert(u);
                            vertex_degree_map->at(u) = compute_core_degree(G, vertex_index_map, u, k, h);

                            if (vertex_degree_map->at(u) < k) {
                                sub_vertex_set->insert(u);
                            }
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }
        auto removed_set = make_shared<unordered_set<uint32_t>>();
        {
            while (!vertex_set->empty()) {
                auto location_vector = pool->split_task(vertex_set);
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);
                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &u = *iter;
                            for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                if (vertex_index_map->at(v)->count(k, h) &&(!vertex_index_map->at(v)->count(k + 1, h +1) || previous_removed_set->count(v))
                                && e_set->size() >= h && !vertex_set->count(v) && !removed_set->count(v)) {
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
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_maintenance::middle_removal_core(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->count(k, h) && (!vertex_index_map->at(u)->count(k + 1, h + 1) || previous_removed_set->count(u))) {
                            sub_clear_set->insert(u);
                            vertex_degree_map->at(u) = compute_core_degree(vertex_edge_index_map, vertex_index_map, u, k, h);

                            if (vertex_degree_map->at(u) < k) {
                                sub_vertex_set->insert(u);
                            }
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }
        auto removed_set = make_shared<unordered_set<uint32_t>>();
        {
            while (!vertex_set->empty()) {
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                auto location_vector = pool->split_task(vertex_set);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                        auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                            auto &u = *iter1;

                            auto edge_map = vertex_edge_index_map->at(u);
                            for (auto iter2 = edge_map->lower_bound(h); iter2 != edge_map->end(); ++iter2) {
                                for (const auto &v: *iter2->second) {
                                    if (vertex_index_map->at(v)->count(k, h) && (!vertex_index_map->at(v)->count(k + 1, h + 1) || previous_removed_set->count(v))
                                        && !vertex_set->count(v) && !removed_set->count(v)) {

                                        vertex_mutex_map->at(v)->lock();
                                        if (vertex_degree_map->at(v) == 0) {
                                            sub_clear_set->insert(v);
                                            vertex_degree_map->at(v) = compute_core_degree(vertex_edge_index_map,
                                                                                           vertex_index_map, v, k, h);
                                        }
                                        --vertex_degree_map->at(v);
                                        vertex_mutex_map->at(v)->unlock();

                                        if (vertex_degree_map->at(v) < k) {
                                            sub_next_vertex_set->insert(v);
                                        }
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
        }

        assign(clear_set, vertex_degree_map, 0, pool);

        return removed_set;
    }


    shared_ptr<unordered_set<uint32_t>>
    branch_multiple_core_maintenance::right_removal_core(const shared_ptr<temporal_graph> &G,
                                                         const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                         const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                         uint32_t k,
                                                         uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &u: *affected_set) {
            if (vertex_index_map->at(u)->get_h(k) == h || previous_removed_set->count(u)) {
                vertex_degree_map->insert({u, compute_core_degree(G, vertex_index_map, u, k, h)});

                if (vertex_degree_map->at(u) < k) {
                    vertex_set->insert(u);
                }
            }
        }
        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if ((vertex_index_map->at(v)->get_h(k) == h || previous_removed_set->count(v))
                        && e_set->size() >= h && !vertex_set->count(v) && !removed_set->count(v)) {
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
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_maintenance::right_removal_core(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
            uint32_t k,
            uint32_t h) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        for (const auto &u: *affected_set) {
            if (vertex_index_map->at(u)->equal(k, h) || previous_removed_set->count(u)) {
                vertex_degree_map->insert({u, compute_core_degree(vertex_edge_index_map, vertex_index_map, u, k, h)});

                if (vertex_degree_map->at(u) < k) {
                    vertex_set->insert(u);
                }
            }
        }
        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!vertex_set->empty()) {
            auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *vertex_set) {
                auto edge_map = vertex_edge_index_map->at(u);
                for (auto iter = edge_map->lower_bound(h); iter != edge_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if ((vertex_index_map->at(v)->equal(k, h) || previous_removed_set->count(v))
                            && !vertex_set->count(v) && !removed_set->count(v)) {
                            if (!vertex_degree_map->count(v)) {
                                vertex_degree_map->insert({v, compute_core_degree(vertex_edge_index_map, vertex_index_map, v, k, h)});
                            }

                            --vertex_degree_map->at(v);
                            if (vertex_degree_map->at(v) < k) {
                                next_vertex_set->insert(v);
                            }
                        }
                    }
                }
            }
            removed_set->merge(*vertex_set);
            swap(*vertex_set, *next_vertex_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>>
    branch_multiple_core_maintenance::right_removal_core(const shared_ptr<temporal_graph> &G,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
                                                         const shared_ptr<unordered_set<uint32_t>> &affected_set,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                         const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
                                                         uint32_t k,
                                                         uint32_t h,
                                                         const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->equal(k, h) || previous_removed_set->count(u)) {
                            sub_clear_set->insert(u);
                            vertex_degree_map->at(u) = compute_core_degree(G, vertex_index_map, u, k, h);

                            if (vertex_degree_map->at(u) < k) {
                                sub_vertex_set->insert(u);
                            }
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        {
            while (!vertex_set->empty()) {
                auto location_vector = pool->split_task(vertex_set);
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                        auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);
                        for (auto iter = sub_begin; iter != sub_end; ++iter) {
                            auto &u = *iter;
                            for (const auto &[v, e_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                if ((vertex_index_map->at(v)->equal(k, h) || previous_removed_set->count(v))
                                    && e_set->size() >= h && !vertex_set->count(v) && !removed_set->count(v)) {
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
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_maintenance::right_removal_core(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &vertex_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &affected_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
            const shared_ptr<unordered_set<uint32_t>> &previous_removed_set,
            uint32_t k,
            uint32_t h,
            const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto clear_set = make_shared<unordered_set<uint32_t>>();
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        if (vertex_index_map->at(u)->equal(k, h) || previous_removed_set->count(u)) {
                            sub_clear_set->insert(u);
                            vertex_degree_map->at(u) = compute_core_degree(vertex_edge_index_map, vertex_index_map, u, k, h);

                            if (vertex_degree_map->at(u) < k) {
                                sub_vertex_set->insert(u);
                            }
                        }
                    }

                    global_mutex->lock();
                    vertex_set->merge(*sub_vertex_set);
                    clear_set->merge(*sub_clear_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }
        auto removed_set = make_shared<unordered_set<uint32_t>>();
        {
            while (!vertex_set->empty()) {
                auto next_vertex_set = make_shared<unordered_set<uint32_t>>();
                auto location_vector = pool->split_task(vertex_set);
                for (uint32_t i = 0; i < thread_number; ++i) {
                    pool->submit_task([=] {
                        auto sub_next_vertex_set = make_shared<unordered_set<uint32_t>>();
                        auto sub_clear_set = make_shared<unordered_set<uint32_t>>();

                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                            auto &u = *iter1;

                            auto edge_map = vertex_edge_index_map->at(u);
                            for (auto iter2 = edge_map->lower_bound(h); iter2 != edge_map->end(); ++iter2) {
                                for (const auto &v: *iter2->second) {
                                    if ((vertex_index_map->at(v)->equal(k, h) || previous_removed_set->count(v))
                                        && !vertex_set->count(v) && !removed_set->count(v)) {

                                        vertex_mutex_map->at(v)->lock();
                                        if (vertex_degree_map->at(v) == 0) {
                                            sub_clear_set->insert(v);
                                            vertex_degree_map->at(v) = compute_core_degree(vertex_edge_index_map,
                                                                                           vertex_index_map, v, k, h);
                                        }
                                        --vertex_degree_map->at(v);
                                        vertex_mutex_map->at(v)->unlock();

                                        if (vertex_degree_map->at(v) < k) {
                                            sub_next_vertex_set->insert(v);
                                        }
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
        }

        assign(clear_set, vertex_degree_map, 0, pool);
        return removed_set;
    }

    void
    branch_multiple_core_maintenance::update_vertex_edge_index_for_insertion(const shared_ptr<scnu::temporal_graph> &G,
                                                                             const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                                             const shared_ptr<unordered_set<uint32_t>>& affected_set) {

        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            affected_set->insert(u);
            affected_set->insert(v);
        }

        for(const auto &u:*affected_set){
            auto u_vertex = G->get_vertex(u);
            if(u_vertex){
                for(const auto &[v, v_edge_set]:*u_vertex->get_neighbor_map()){
                    if(affected_set->count(v)){
                        auto h = v_edge_set->size();
                        vertex_edge_index_map->at(u)->at(h)->erase(v);
                        if(vertex_edge_index_map->at(u)->at(h)->empty()){
                            vertex_edge_index_map->at(u)->erase(h);
                        }
                    }
                }
            }else{
                vertex_edge_index_map->insert({u, make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
            }
        }

        G->insert_edge_collection(edge_set);

        for(const auto &u:*affected_set){
            for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                if(affected_set->count(v)){
                    auto h = v_edge_set->size();
                    if (!vertex_edge_index_map->at(u)->count(h)) {
                        vertex_edge_index_map->at(u)->insert({h, make_shared<unordered_set<uint32_t>>()});
                    }
                    vertex_edge_index_map->at(u)->at(h)->insert(v);
                }
            }
        }
    }

    void
    branch_multiple_core_maintenance::update_vertex_edge_index_for_insertion(const shared_ptr<scnu::temporal_graph> &G,
                                                                             const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                                             const shared_ptr<unordered_set<uint32_t>>& affected_set,
                                                                             const shared_ptr<thread_pool>& pool) {

        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            affected_set->insert(u);
            affected_set->insert(v);
        }

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();
        {
            auto location_vector = pool->split_task(affected_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto new_vertex_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &u = *iter;

                        auto u_vertex = G->get_vertex(u);
                        if(u_vertex){
                            for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                                if(affected_set->count(v)){
                                    auto h = v_edge_set->size();
                                    vertex_edge_index_map->at(u)->at(h)->erase(v);
                                    if(vertex_edge_index_map->at(u)->at(h)->empty()){
                                        vertex_edge_index_map->at(u)->erase(h);
                                    }
                                }
                            }
                        }else{
                            new_vertex_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    for(const auto &u:*new_vertex_set){
                        vertex_edge_index_map->insert({u, make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
                    }
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        G->insert_edge_collection(edge_set);
        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                            if(affected_set->count(v)){
                                auto h = v_edge_set->size();
                                if (!vertex_edge_index_map->at(u)->count(h)) {
                                    vertex_edge_index_map->at(u)->insert({h, make_shared<unordered_set<uint32_t>>()});
                                }
                                vertex_edge_index_map->at(u)->at(h)->insert(v);
                            }
                        }
                    }
                });
            }
            pool->barrier();
        }
    }

    void
    branch_multiple_core_maintenance::update_vertex_edge_index_for_removal(const shared_ptr<scnu::temporal_graph> &G,
                                                                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                                           const shared_ptr<unordered_set<uint32_t>>& affected_set,
                                                                           const shared_ptr<unordered_set<uint32_t>>& isolated_vertex_set) {
        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            affected_set->insert(u);
            affected_set->insert(v);
        }

        for(const auto &u:*affected_set){
            for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                if(affected_set->count(v)){
                    auto h = v_edge_set->size();
                    vertex_edge_index_map->at(u)->at(h)->erase(v);
                    if(vertex_edge_index_map->at(u)->at(h)->empty()){
                        vertex_edge_index_map->at(u)->erase(h);
                    }
                }
            }
        }

        G->remove_edge_collection(edge_set, isolated_vertex_set);

        for(const auto &u:*isolated_vertex_set){
            affected_set->erase(u);
        }

        for(const auto &u:*affected_set){
            for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                if(affected_set->count(v)){
                    auto h = v_edge_set->size();
                    if (!vertex_edge_index_map->at(u)->count(h)) {
                        vertex_edge_index_map->at(u)->insert({h, make_shared<unordered_set<uint32_t>>()});
                    }
                    vertex_edge_index_map->at(u)->at(h)->insert(v);
                }
            }
        }
    }

    void
    branch_multiple_core_maintenance::update_vertex_edge_index_for_removal(const shared_ptr<scnu::temporal_graph> &G,
                                                                           const shared_ptr<unordered_set<shared_ptr<temporal_edge>>> &edge_set,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &vertex_edge_index_map,
                                                                           const shared_ptr<unordered_set<uint32_t>>& affected_set,
                                                                           const shared_ptr<unordered_set<uint32_t>>& isolated_vertex_set,
                                                                           const shared_ptr<thread_pool>& pool) {
        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            affected_set->insert(u);
            affected_set->insert(v);
        }

        auto thread_number = pool->get_thread_number();
        {
            auto location_vector = pool->split_task(affected_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter != sub_end; ++iter){
                        auto &u = *iter;

                        for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                            if(affected_set->count(v)){
                                auto h = v_edge_set->size();
                                vertex_edge_index_map->at(u)->at(h)->erase(v);
                                if(vertex_edge_index_map->at(u)->at(h)->empty()){
                                    vertex_edge_index_map->at(u)->erase(h);
                                }
                            }
                        }
                    }
                });
            }
            pool->barrier();
        }

        G->remove_edge_collection(edge_set, isolated_vertex_set);

        for(const auto &u:*isolated_vertex_set){
            affected_set->erase(u);
        }

        {
            auto location_vector = pool->split_task(affected_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;
                        for(const auto &[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                            if(affected_set->count(v)){
                                auto h = v_edge_set->size();
                                if (!vertex_edge_index_map->at(u)->count(h)) {
                                    vertex_edge_index_map->at(u)->insert({h, make_shared<unordered_set<uint32_t>>()});
                                }
                                vertex_edge_index_map->at(u)->at(h)->insert(v);
                            }
                        }
                    }
                });
            }
            pool->barrier();
        }
    }

}