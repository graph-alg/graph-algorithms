
#include "multiple_core/basic_multiple_core_decomposition.h"

namespace scnu {

    void basic_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
        }
    }


    void basic_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i< thread_number;++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i+1);
                for(auto iter = sub_begin; iter!=sub_end;++iter){
                    auto u = iter->first;
                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                }
            });
        }
        pool->barrier();
    }


    void basic_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
            vertex_mutex_map->insert({u, shared_ptr<mutex>()});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i< thread_number;++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i+1);
                for(auto iter = sub_begin; iter!=sub_end;++iter){
                    auto &u = iter->first;
                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                    vertex_mutex_map->at(u) = make_shared<mutex>();
                }
            });
        }
        pool->barrier();
    }

    void basic_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph> &G,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->at(u)->insert(1, 1);
            vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
        }

        for (uint32_t h = 1; !vertex_degree_map->empty(); ++h) {
            auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);

            for (uint32_t k = 1; !sub_vertex_degree_map->empty(); ++k) {
                auto removed_set = find_k_core(G, sub_vertex_degree_map, k + 1, h);
                for (const auto &u: *removed_set) {
                    for (uint64_t index = 1; index <= k; ++index) {
                        vertex_index_map->at(u)->insert(index, h);
                    }
                }
            }

            auto removed_set = find_h_core(G, vertex_degree_map, h + 1);
            for (const auto &u: *removed_set) {
                vertex_index_map->at(u)->insert(1, h);
            }
        }
    }

    void basic_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph> &G,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                      const shared_ptr<thread_pool>& pool) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
        }

        auto global_mutex = make_shared<mutex>();
        for (uint32_t h = 1; !vertex_degree_map->empty(); ++h) {
            auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);

            pool->submit_task([=] {
                for (uint32_t k = 1; !sub_vertex_degree_map->empty(); ++k) {
                    auto removed_set = find_k_core(G, sub_vertex_degree_map, k + 1, h);

                    global_mutex->lock();
                    for (const auto &u: *removed_set) {
                        for (uint64_t index = 1; index <= k; ++index) {
                            vertex_index_map->at(u)->insert(index, h);
                        }
                    }
                    global_mutex->unlock();
                }
            });

            auto removed_set = find_h_core(G, vertex_degree_map, h + 1);
            global_mutex->lock();
            for (const auto &u: *removed_set) {
                vertex_index_map->at(u)->insert(1, h);
            }
            global_mutex->unlock();
        }
        pool->barrier();
    }

    void basic_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph> &G,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> & vertex_mutex_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                      const shared_ptr<thread_pool>& pool) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
        }

        for (uint32_t h = 1; !vertex_degree_map->empty(); ++h) {
            auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
            for (uint32_t k = 1; !sub_vertex_degree_map->empty(); ++k) {
                auto removed_set = find_k_core(G, vertex_mutex_map, sub_vertex_degree_map, k + 1, h, pool);
                assign(removed_set, vertex_index_map, k, h,pool);
            }
            auto removed_set = find_h_core(G, vertex_degree_map, h + 1, pool);
            assign(removed_set,vertex_index_map, 1, h, pool);
        }
    }

    void basic_multiple_core_decomposition::assign(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                   uint32_t k,
                                                   uint32_t h,
                                                   const shared_ptr<thread_pool>& pool){
        if (vertex_set->empty()) {
            return;
        }
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_set);
        for(uint32_t i = 0; i < thread_number;++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto u = *iter;
                    for (uint64_t index = 1; index <= k; ++index) {
                        vertex_index_map->at(u)->insert(index, h);
                    }
                }
            });
        }
        pool->barrier();
    }


    shared_ptr<unordered_set<uint32_t>> basic_multiple_core_decomposition::find_h_core(const shared_ptr<temporal_graph> &G,
                                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                                       uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for (auto iter = vertex_degree_map->begin(); iter != vertex_degree_map->end();) {
            auto &[u, u_degree] =*iter;
            ++iter;
            for (const auto&[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                if (vertex_degree_map->count(v) && v_edge_set->size() == h - 1) {
                    --vertex_degree_map->at(u);
                }
            }
            if (vertex_degree_map->at(u) < 1) {
                evicted_set->insert(u);
            }
        }

        for(const auto &u:*evicted_set){
            vertex_degree_map->erase(u);
        }

        return evicted_set;
    }

    shared_ptr<unordered_set<uint32_t>> basic_multiple_core_decomposition::find_h_core(const shared_ptr<temporal_graph> &G,
                                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                                       uint32_t h,
                                                                                       const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();

        auto global_mutex = make_shared<mutex>();
        auto location_vector = pool->split_task(vertex_degree_map);
        for(uint32_t i = 0; i < thread_number;++i){
            pool->submit_task([=]{
                auto sub_evicted_set =  make_shared<unordered_set<uint32_t>>();

                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                    auto &[u, u_degree] = *iter;
                    for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                        if (vertex_degree_map->count(v) && v_edge_set->size() == h - 1) {
                            --vertex_degree_map->at(u);
                        }
                    }
                    if (vertex_degree_map->at(u) < 1) {
                        sub_evicted_set->insert(u);
                    }
                }

                global_mutex->lock();
                evicted_set->merge(*sub_evicted_set);
                global_mutex->unlock();
            });
        }
        pool->barrier();
        for(const auto &u:*evicted_set){
            vertex_degree_map->erase(u);
        }
        return evicted_set;
    }


    shared_ptr<unordered_set<uint32_t>> basic_multiple_core_decomposition::find_k_core(const shared_ptr<temporal_graph> &G,
                                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                                       uint32_t k,
                                                                                       uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for (const auto&[u, u_degree]: *vertex_degree_map) {
            if (u_degree < k) {
                evicted_set->insert(u);
            }
        }

        if(evicted_set->size() == vertex_degree_map->size()){
            vertex_degree_map->clear();
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
            for (const auto &u: *evicted_set) {
                vertex_degree_map->erase(u);
                removed_set->insert(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *evicted_set) {
                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_degree_map->count(v) && v_edge_set->size() >= h && vertex_degree_map->at(v) >= k ) {
                        --vertex_degree_map->at(v);
                        if (vertex_degree_map->at(v) < k) {
                            next_evicted_set->insert(v);
                        }
                    }
                }
            }
            swap(*evicted_set, *next_evicted_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> basic_multiple_core_decomposition::find_k_core(const shared_ptr<temporal_graph> &G,
                                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> & vertex_mutex_map,
                                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                                                       uint32_t k,
                                                                                       uint32_t h,
                                                                                       const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto global_mutex = make_shared<mutex>();
        {
            auto location_vector = pool->split_task(vertex_degree_map);
            for(uint32_t i = 0; i < thread_number;++i){
                pool->submit_task([=]{
                    auto sub_vertex_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &[u, u_degree] = *iter;
                        if (vertex_degree_map->at(u) < k) {
                            sub_vertex_set->insert(u);
                        }
                    }

                    global_mutex->lock();
                    evicted_set->merge(*sub_vertex_set);
                    global_mutex->unlock();
                });
            }
            pool->barrier();
        }

        if(evicted_set->size() == vertex_degree_map->size()){
            vertex_degree_map->clear();
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
            for(const auto&u:*evicted_set){
                vertex_degree_map->erase(u);
                removed_set->insert(u);
            }
            auto location_vector = pool->split_task(evicted_set);
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for(uint32_t i = 0; i < thread_number;++i){
                pool->submit_task([=]{
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_removed_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;

                        for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (vertex_degree_map->count(v) && v_edge_set->size() >= h && vertex_degree_map->at(v) >= k){

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
            swap(*evicted_set, *next_evicted_set);
        }
        return removed_set;
    }
}
