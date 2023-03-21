
#include "multiple_core/hierarchy_multiple_core_decomposition.h"

namespace scnu{

    void hierarchy_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->insert({u, make_shared<multiple_core_pair_map_index>()});
        }
    }


    void hierarchy_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                 const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i < thread_number;++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter!=sub_end;++iter){
                    auto &u = iter->first;
                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                }
            });
        }
        pool->barrier();
    }


    void hierarchy_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                     const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_mutex_map->insert({u, shared_ptr<mutex>()});
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter != sub_end; ++iter){
                    auto &u = iter->first;
                    vertex_mutex_map->at(u) = make_shared<mutex>();
                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                }
            });
        }
        pool->barrier();
    }

    void hierarchy_multiple_core_decomposition::assign(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                        uint32_t k,
                                                        uint32_t h,
                                                        const shared_ptr<thread_pool>& pool){
        if(vertex_set->empty()){
            return;
        }
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for (auto iter = sub_begin; iter!=sub_end;++iter) {
                    auto &u = *iter;
                    for (uint64_t index = 1; index <= k; ++index) {
                        vertex_index_map->at(u)->insert(index, h);
                    }
                }
            });
        }
        pool->barrier();
    }


    void hierarchy_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph>& G,
                                                          const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map)
    {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->at(u)->insert(1, 1);
            vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
        }

        auto current_layer_vertex_degree_map = make_shared<unordered_map<pair<uint32_t,uint32_t>,
                shared_ptr<unordered_map<uint32_t, uint32_t>>,hash_pair,equal_pair>>();
        current_layer_vertex_degree_map->insert({{1, 1}, vertex_degree_map});


        while (!current_layer_vertex_degree_map->empty()) {
            auto next_layer_vertex_degree_map = make_shared<unordered_map<pair<uint32_t,uint32_t>,
                    shared_ptr<unordered_map<uint32_t, uint32_t>>,hash_pair,equal_pair>>();
            for(const auto&p:*current_layer_vertex_degree_map)
            {
                auto [k, h] = p.first;
                auto current_vertex_degree_map = p.second;

                if(k == 1){
                    auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*current_vertex_degree_map);
                    auto removed_set = find_h_core(G, sub_vertex_degree_map, h + 1);
                    for (const auto& u: *removed_set) {
                        vertex_index_map->at(u)->insert(1, h);
                    }
                    if(!sub_vertex_degree_map->empty()){
                        next_layer_vertex_degree_map->insert({{1, h + 1}, sub_vertex_degree_map});
                    }
                }

                auto &sub_vertex_degree_map = current_vertex_degree_map;
                auto removed_set = find_k_core(G, sub_vertex_degree_map, k + 1, h);
                for (const auto&u: *removed_set) {
                    for (uint64_t index = 1; index <= k; ++index) {
                        vertex_index_map->at(u)->insert(index, h);
                    }
                }
                if(!sub_vertex_degree_map->empty()){
                    next_layer_vertex_degree_map->insert({{k + 1, h}, sub_vertex_degree_map});
                }
            }
            swap(*current_layer_vertex_degree_map, *next_layer_vertex_degree_map);
        }
    }

    void hierarchy_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph>& G,
                                                          const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                          const shared_ptr<thread_pool>& pool){
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
        }

        auto current_layer_vertex_degree_map = make_shared<unordered_map<pair<uint32_t,uint32_t>,
                shared_ptr<unordered_map<uint32_t, uint32_t>>,hash_pair,equal_pair>>();
        current_layer_vertex_degree_map->insert({{1, 1}, vertex_degree_map});

        auto global_mutex = make_shared<mutex>();
        while (!current_layer_vertex_degree_map->empty()) {
            auto next_layer_vertex_degree_map = make_shared<unordered_map<pair<uint32_t,uint32_t>,
                    shared_ptr<unordered_map<uint32_t, uint32_t>>,hash_pair,equal_pair>>();
            if(current_layer_vertex_degree_map->size() > 1){
                for (const auto &p: *current_layer_vertex_degree_map) {
                    auto [k, h] = p.first;
                    auto current_vertex_degree_map = p.second;
                    if (k == 1) {
                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*current_vertex_degree_map);
                        {
                            pool->submit_task([=] {
                                auto removed_set = find_h_core(G, sub_vertex_degree_map, h + 1);

                                global_mutex->lock();
                                for (const auto &u: *removed_set) {
                                    vertex_index_map->at(u)->insert(1, h);
                                }
                                if(!sub_vertex_degree_map->empty()){
                                    next_layer_vertex_degree_map->insert({{1, h + 1}, sub_vertex_degree_map});
                                }
                                global_mutex->unlock();
                            });
                        }
                    }
                    {
                        auto &sub_vertex_degree_map = current_vertex_degree_map;
                        pool->submit_task([=] {
                            auto removed_set = find_k_core(G, sub_vertex_degree_map, k + 1, h);

                            global_mutex->lock();
                            for (const auto& u: *removed_set) {
                                for (uint64_t index = 1; index <= k; ++index) {
                                    vertex_index_map->at(u)->insert(index, h);
                                }
                            }
                            if(!sub_vertex_degree_map->empty()){
                                next_layer_vertex_degree_map->insert({{k + 1, h}, sub_vertex_degree_map});
                            }
                            global_mutex->unlock();
                        });
                    }
                }
                pool->barrier();
            }else{
                for(const auto&p:*current_layer_vertex_degree_map)
                {
                    auto [k, h] = p.first;
                    auto current_vertex_degree_map = p.second;

                    if(k == 1){
                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*current_vertex_degree_map);
                        auto removed_set = find_h_core(G, sub_vertex_degree_map, h + 1);
                        for (const auto& u: *removed_set) {
                            vertex_index_map->at(u)->insert(1, h);
                        }
                        if(!sub_vertex_degree_map->empty()){
                            next_layer_vertex_degree_map->insert({{1, h + 1}, sub_vertex_degree_map});
                        }
                    }

                    auto &sub_vertex_degree_map = current_vertex_degree_map;
                    auto removed_set = find_k_core(G, sub_vertex_degree_map, k + 1, h);
                    for (const auto&u: *removed_set) {
                        for (uint64_t index = 1; index <= k; ++index) {
                            vertex_index_map->at(u)->insert(index, h);
                        }
                    }
                    if(!sub_vertex_degree_map->empty()){
                        next_layer_vertex_degree_map->insert({{k + 1, h}, sub_vertex_degree_map});
                    }
                }
            }

            swap(*current_layer_vertex_degree_map, *next_layer_vertex_degree_map);
        }
    }

    void hierarchy_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph>& G,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> & vertex_mutex_map,
                                                          const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                          const shared_ptr<thread_pool>& pool){
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->at(u)->insert(1, 1);
            vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
        }

        auto current_layer_vertex_degree_map = make_shared<unordered_map<pair<uint32_t,uint32_t>,
                shared_ptr<unordered_map<uint32_t, uint32_t>>,hash_pair,equal_pair>>();
        current_layer_vertex_degree_map->insert({{1, 1}, vertex_degree_map});

        while (!current_layer_vertex_degree_map->empty()) {
            auto next_layer_vertex_degree_map = make_shared<unordered_map<pair<uint32_t,uint32_t>,
                    shared_ptr<unordered_map<uint32_t, uint32_t>>,hash_pair,equal_pair>>();
            for (const auto &p: *current_layer_vertex_degree_map) {
                auto[k, h] = p.first;
                auto current_vertex_degree_map = p.second;

                if (k == 1){
                    auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*current_vertex_degree_map);

                    auto removed_set = find_h_core(G, vertex_mutex_map, sub_vertex_degree_map, h + 1, pool);
                    assign(removed_set,vertex_index_map,k, h,pool);
                    if (!sub_vertex_degree_map->empty()) {
                        next_layer_vertex_degree_map->insert({{1, h + 1}, sub_vertex_degree_map});
                    }
                }
                {
                    auto &sub_vertex_degree_map = current_vertex_degree_map;
                    auto removed_set = find_k_core(G, vertex_mutex_map, sub_vertex_degree_map, k + 1, h, pool);
                    assign(removed_set, vertex_index_map,k, h,pool);
                    if (!sub_vertex_degree_map->empty()) {
                        next_layer_vertex_degree_map->insert({{k + 1, h}, sub_vertex_degree_map});
                    }
                }
            }
            swap(*current_layer_vertex_degree_map, *next_layer_vertex_degree_map);
        }
    }



    shared_ptr<unordered_set<uint32_t>> hierarchy_multiple_core_decomposition::find_h_core(const shared_ptr<temporal_graph> &G,
                                                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                           uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for (const auto&[u, u_degree]: *vertex_degree_map) {
            for (const auto&[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                if (vertex_degree_map->count(v) && v_edge_set->size() == h - 1) {
                    --vertex_degree_map->at(u);
                }
            }
            if (vertex_degree_map->at(u) < 1) {
                evicted_set->insert(u);
            }
        }

        for (const auto &v: *evicted_set) {
            vertex_degree_map->erase(v);
        }
        return evicted_set;
    }

    shared_ptr<unordered_set<uint32_t>> hierarchy_multiple_core_decomposition::find_h_core(const shared_ptr<temporal_graph> &G,
                                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> & vertex_mutex_map,
                                                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                           uint32_t h,
                                                                                           const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto global_mutex = make_shared<mutex>();
        {
            auto location_vector = pool->split_task(vertex_degree_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
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
        }

        for (const auto &v: *evicted_set) {
            vertex_degree_map->erase(v);
        }
        return evicted_set;
    }

    shared_ptr<unordered_set<uint32_t>> hierarchy_multiple_core_decomposition::find_k_core(const shared_ptr<temporal_graph> &G,
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
            for(const auto &u:*evicted_set) {
                vertex_degree_map->erase(u);
                removed_set->insert(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*evicted_set){
                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_degree_map->count(v) && v_edge_set->size() >= h && vertex_degree_map->at(v) >= k) {
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

    shared_ptr<unordered_set<uint32_t>> hierarchy_multiple_core_decomposition::find_k_core(const shared_ptr<temporal_graph> &G,
                                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> & vertex_mutex_map,
                                                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                           uint32_t k,
                                                                                           uint32_t h,
                                                                                           const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        auto global_mutex = make_shared<mutex>();
        {
            auto location_vector = pool->split_task(vertex_degree_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto sub_begin = *location_vector->at(i);
                    auto sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &[u, u_degree] = *iter;
                        if (u_degree < k) {
                            sub_evicted_set->insert(u);
                        }
                    }
                    global_mutex->lock();
                    evicted_set->merge(*sub_evicted_set);
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
            for(const auto &u:*evicted_set){
                vertex_degree_map->erase(u);
                removed_set->insert(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(evicted_set);
            for(uint32_t i = 0; i < thread_number;++i){
                pool->submit_task([=]{
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end;  ++iter) {
                        auto &u = *iter;
                        for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (vertex_degree_map->count(v) && v_edge_set->size() >= h && vertex_degree_map->at(v) >= k) {

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