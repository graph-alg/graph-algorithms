
#include "multiple_core/branch_multiple_core_decomposition.h"

namespace scnu {

    void branch_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter != sub_end; ++iter){
                    auto u = iter->first;
                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                }
            });
        }
        pool->barrier();
    }

    void branch_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t,
                                                          shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
            vertex_edge_size_map->insert({u, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter != sub_end; ++iter){
                    auto u = iter->first;

                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                    vertex_edge_size_map->at(u) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();

                    auto u_vertex = G->get_vertex(u);
                    for (const auto&[v, v_edge_set]: *u_vertex->get_neighbor_map()) {
                        auto h = v_edge_set->size();
                        if(!vertex_edge_size_map->at(u)->count(h)){
                            vertex_edge_size_map->at(u)->insert({h, make_shared<unordered_set<uint32_t>>()});
                        }
                        vertex_edge_size_map->at(u)->at(h)->insert(v);
                    }
                }
            });
        }
        pool->barrier();
    }

    void branch_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_mutex_map->insert({u, make_shared<mutex>()});
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
        }

        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter != sub_end && iter!=vertex_index_map->end(); ++iter){
                    auto u = iter->first;
                    vertex_mutex_map->at(u) = make_shared<mutex>();
                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                }
            });
        }
        pool->barrier();
    }

    void branch_multiple_core_decomposition::init(const shared_ptr<scnu::temporal_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                  const shared_ptr<unordered_map<uint32_t,
                                                          shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::multiple_core_pair_map_index>>> &vertex_index_map,
                                                  const shared_ptr<thread_pool>& pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            vertex_mutex_map->insert({u, shared_ptr<mutex>()});
            vertex_index_map->insert({u, shared_ptr<multiple_core_pair_map_index>()});
            vertex_edge_size_map->insert({u, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_index_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);
                for(auto iter = sub_begin; iter != sub_end && iter!=vertex_index_map->end(); ++iter){
                    auto u = iter->first;

                    vertex_mutex_map->at(u) = make_shared<mutex>();
                    vertex_index_map->at(u) = make_shared<multiple_core_pair_map_index>();
                    vertex_edge_size_map->at(u) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();

                    auto u_vertex = G->get_vertex(u);
                    for (const auto&[v, v_edge_set]: *u_vertex->get_neighbor_map()) {
                        auto h = v_edge_set->size();
                        if(!vertex_edge_size_map->at(u)->count(h)){
                            vertex_edge_size_map->at(u)->insert({h, make_shared<unordered_set<uint32_t>>()});
                        }
                        vertex_edge_size_map->at(u)->at(h)->insert(v);
                    }
                }
            });
        }
        pool->barrier();
    }

    uint32_t branch_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph>& G,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                           const shared_ptr<thread_pool>& pool){
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        {
            for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
                vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
            }
        }

        uint32_t max_delta = 0;

        auto global_mutex = make_shared<mutex>();
        for(uint32_t delta = 1; !vertex_degree_map->empty(); ++delta) {
            max_delta = delta;
            auto current_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
            pool->submit_task([=] {
                if(delta == 1){
                    auto &sub_vertex_degree_map = current_vertex_degree_map;
                    /**
                     * @brief compute (1, h)-cores
                     */
                    auto h_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                    for (const auto &[u, u_degree]: *sub_vertex_degree_map) {
                        uint32_t max_h = 1;
                        for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (v_edge_set->size() > max_h) {
                                max_h = v_edge_set->size();
                            }
                        }
                        if(!h_map->count(max_h)){
                            h_map->insert({max_h, make_shared<unordered_set<uint32_t>>()});
                        }
                        h_map->at(max_h)->insert(u);
                    }
                    global_mutex->lock();
                    for(const auto &[h, u_set]:*h_map){
                        for(const auto&u:*u_set){
                            vertex_index_map->at(u)->insert(1, h);
                        }
                    }
                    global_mutex->unlock();
                    /**
                     * @brief compute (k,1)-cores
                     */
                    for (uint32_t k = 1; !sub_vertex_degree_map->empty(); ++k) {
                        auto removed_set = find_k_core(G, sub_vertex_degree_map, k + 1, 1);
                        global_mutex->lock();
                        for (const auto &u: *removed_set) {
                            vertex_index_map->at(u)->insert(k, 1);
                        }
                        global_mutex->unlock();
                    }
                }else {
                    {
                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*current_vertex_degree_map);
                        for (uint32_t k = delta + 1; !sub_vertex_degree_map->empty(); ++k) {
                            auto removed_set = find_k_core(G, sub_vertex_degree_map, k, delta);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                for (uint32_t index = delta; index < k; ++index) {
                                    vertex_index_map->at(u)->insert(index, delta);
                                }
                            }
                            global_mutex->unlock();
                        }
                    }

                    {
                        auto &sub_vertex_map = current_vertex_degree_map;
                        for (uint32_t h = delta + 1; !sub_vertex_map->empty(); ++h) {
                            auto removed_set = find_h_core(G, sub_vertex_map, delta, h);

                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                for (uint32_t index = delta; index < h; ++index) {
                                    vertex_index_map->at(u)->insert(delta, index);
                                }
                            }
                            global_mutex->unlock();
                        }
                    }
                }
            });
            {
                auto removed_set = find_core(G, vertex_degree_map, delta + 1, delta + 1);
                global_mutex->lock();
                for (const auto &u: *removed_set) {
                    vertex_index_map->at(u)->insert(delta, delta);
                }
                global_mutex->unlock();
            }
        }
        pool->barrier();
        return max_delta;
    }

    uint32_t branch_multiple_core_decomposition:: decompose(const shared_ptr<temporal_graph>& G,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                            const shared_ptr<thread_pool>& pool)
    {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        {
            for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
                vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
            }
        }

        uint32_t max_delta = 0;
        auto global_mutex = make_shared<mutex>();
        for(uint32_t delta = 1; !vertex_degree_map->empty(); ++delta) {
            max_delta = delta;
            auto current_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
            pool->submit_task([=] {
                if(delta == 1){
                    auto &sub_vertex_degree_map = current_vertex_degree_map;
                    /**
                     * @brief compute (1, h)-cores
                     */
                    auto h_map = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                    for (const auto &[u, u_degree]: *sub_vertex_degree_map) {
                        auto  max_h = vertex_edge_index_map->at(u)->rbegin()->first;
                        if(!h_map->count(max_h)){
                            h_map->insert({max_h, make_shared<unordered_set<uint32_t>>()});
                        }
                        h_map->at(max_h)->insert(u);
                    }
                    global_mutex->lock();
                    for(const auto &[h, u_set]:*h_map){
                        for(const auto&u:*u_set){
                            vertex_index_map->at(u)->insert(1, h);
                        }
                    }
                    global_mutex->unlock();
                    /**
                     * @brief compute (k,1)-cores
                     */
                    for (uint32_t k = 1; !sub_vertex_degree_map->empty(); ++k) {
                        auto removed_set = find_k_core(vertex_edge_index_map, sub_vertex_degree_map, k + 1, 1);
                        global_mutex->lock();
                        for (const auto &u: *removed_set) {
                            vertex_index_map->at(u)->insert(k, 1);
                        }
                        global_mutex->unlock();
                    }
                }else {
                    {
                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*current_vertex_degree_map);
                        for (uint32_t k = delta + 1; !sub_vertex_degree_map->empty(); ++k) {
                            auto removed_set = find_k_core(vertex_edge_index_map, sub_vertex_degree_map, k, delta);
                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                for (uint32_t index = delta; index < k; ++index) {
                                    vertex_index_map->at(u)->insert(index, delta);
                                }
                            }
                            global_mutex->unlock();
                        }
                    }

                    {
                        auto &sub_vertex_degree_map = current_vertex_degree_map;
                        auto h_map = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                        for (uint32_t h = delta + 1; !sub_vertex_degree_map->empty(); ++h) {
                            auto removed_set = find_h_core(vertex_edge_index_map, sub_vertex_degree_map, delta, h);
                            global_mutex->lock();
                            for (const auto &u: *removed_set) {
                                for (uint32_t index = delta; index < h; ++index) {
                                    vertex_index_map->at(u)->insert(delta, index);
                                }
                            }
                            global_mutex->unlock();
                        }
                    }
                }
            });
            {
                auto removed_set = find_core(vertex_edge_index_map, vertex_degree_map, delta + 1, delta + 1);
                global_mutex->lock();
                for (const auto &u: *removed_set) {
                    vertex_index_map->at(u)->insert(delta, delta);
                }
                global_mutex->unlock();
            }
        }
        pool->barrier();
        return max_delta;
    }

    uint32_t branch_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph>& G,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                           const shared_ptr<thread_pool>& pool){
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        {
            for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
                vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
            }
        }

        uint32_t max_delta = 0;

        for(uint32_t delta = 1; !vertex_degree_map->empty(); ++delta) {
            max_delta = delta;
            if (delta == 1) {
                /**
                 * @brief compute (1, h)-cores
                 */
                auto location_vector = pool->split_task(vertex_degree_map);
                auto thread_number = pool->get_thread_number();
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for(auto iter = sub_begin; iter != sub_end; ++iter){
                            auto &u = iter->first;
                            uint32_t  max_h = 1;
                            for(const auto&[v, v_edge_set]:*G->get_vertex(u)->get_neighbor_map()){
                                if(v_edge_set->size() > max_h){
                                    max_h = v_edge_set->size();
                                }
                            }
                            vertex_index_map->at(u)->insert(1, max_h);
                        }
                    });
                }
                auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);;
                pool->barrier();
                /**
                 * @brief compute (k,1)-cores
                 */
                for (uint32_t k = 1; !sub_vertex_degree_map->empty(); ++k) {
                    auto removed_set = find_k_core(G, vertex_mutex_map, sub_vertex_degree_map, k + 1, 1, pool);
                    assign(removed_set, vertex_index_map, k, 1, pool);
                }
            } else {
                auto sub_vertex_degree_map1 = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto sub_vertex_degree_map2 = make_shared<unordered_map<uint32_t, uint32_t>>();
                pool->submit_task([=]{
                    std::copy(vertex_degree_map->begin(), vertex_degree_map->end(), std::inserter(*sub_vertex_degree_map1, sub_vertex_degree_map1->end()));
                });
                std::copy(vertex_degree_map->begin(), vertex_degree_map->end(), std::inserter(*sub_vertex_degree_map2, sub_vertex_degree_map2->end()));
                pool->barrier();
                {
                    auto &sub_vertex_degree_map = sub_vertex_degree_map1;
                    for (uint32_t k = delta + 1; !sub_vertex_degree_map->empty(); ++k) {
                        auto removed_set = find_k_core(G, vertex_mutex_map, sub_vertex_degree_map, k, delta, pool);
                        assign(removed_set, vertex_index_map, k - 1, delta, pool);
                    }
                }

                {
                    auto &sub_vertex_degree_map = sub_vertex_degree_map2;
                    for (uint32_t h = delta + 1; !sub_vertex_degree_map->empty(); ++h) {
                        auto removed_set = find_h_core(G, vertex_mutex_map, sub_vertex_degree_map, delta, h, pool);
                        assign(removed_set, vertex_index_map, delta, h - 1, pool);
                    }
                }
            }
            {
                auto removed_set = find_core(G, vertex_mutex_map,vertex_degree_map, delta + 1, delta + 1, pool);
                assign(removed_set, vertex_index_map, delta, delta, pool);
            }
        }
        return max_delta;
    }

    uint32_t branch_multiple_core_decomposition::decompose(const shared_ptr<temporal_graph>& G,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_index_map,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                           const shared_ptr<thread_pool>& pool)
    {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        {
            for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
                vertex_degree_map->insert({u, u_vertex->get_neighbor_size()});
            }
        }

        uint32_t max_delta = 0;

        for(uint32_t delta = 1; !vertex_degree_map->empty(); ++delta) {
            max_delta = delta;
            if (delta == 1) {

                /**
                 * @brief compute (1, h)-cores
                 */
                auto location_vector = pool->split_task(vertex_degree_map);
                auto thread_number = pool->get_thread_number();
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);

                        for(auto iter = sub_begin; iter != sub_end; ++iter){
                            auto &u = iter->first;
                            auto  max_h = vertex_edge_index_map->at(u)->rbegin()->first;
                            vertex_index_map->at(u)->insert(1, max_h);
                        }
                    });
                }
                auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
                pool->barrier();
                /**
                 * @brief compute (k,1)-cores
                 */
                for (uint32_t k = 1; !sub_vertex_degree_map->empty(); ++k) {
                    auto removed_set = find_k_core(vertex_edge_index_map, vertex_mutex_map, sub_vertex_degree_map, k + 1, 1, pool);
                    assign(removed_set, vertex_index_map, k, 1, pool);
                }
            } else {
                {
                    auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
                    for (uint32_t k = delta + 1; !sub_vertex_degree_map->empty(); ++k) {
                        auto removed_set = find_k_core(vertex_edge_index_map, vertex_mutex_map, sub_vertex_degree_map, k, delta, pool);
                        assign(removed_set, vertex_index_map, k - 1, delta, pool);
                    }
                }

                {
                    auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map);
                    for (uint32_t h = delta + 1; !sub_vertex_degree_map->empty(); ++h) {
                        auto removed_set = find_h_core(vertex_edge_index_map, vertex_mutex_map, sub_vertex_degree_map, delta, h, pool);
                        assign(removed_set, vertex_index_map, delta, h - 1, pool);
                    }
                }
            }
            {
                auto removed_set = find_core(vertex_edge_index_map, vertex_mutex_map,vertex_degree_map, delta + 1, delta + 1, pool);
                assign(removed_set, vertex_index_map, delta, delta, pool);
            }
        }
        return max_delta;
    }

    void branch_multiple_core_decomposition::assign(const shared_ptr<unordered_set<uint32_t>>& vertex_set,
                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<multiple_core_pair_map_index>>>& vertex_index_map,
                                                    uint32_t k,
                                                    uint32_t h,
                                                    const shared_ptr<thread_pool>& pool){
        if (vertex_set->empty()) {
            return;
        }
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(vertex_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_begin = *location_vector->at(i);
                auto sub_end = *location_vector->at(i + 1);

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto u = *iter;
                    if(k > h){
                        for(uint32_t index = h; index <= k;++index){
                            vertex_index_map->at(u)->insert(index, h);
                        }
                    }else if(k < h){
                        for(uint32_t index = k; index <= h;++index){
                            vertex_index_map->at(u)->insert(k, index);
                        }
                    }else if(k == h || h == 1){
                        vertex_index_map->at(u)->insert(k, h);
                    }
                }
            });
        }
        pool->barrier();
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_core(const shared_ptr<temporal_graph>& G,
                                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                                                      uint32_t k,
                                                                                      uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for(const auto &[u, u_degree]:*vertex_map) {
            if(vertex_map->at(u) < k){
                evicted_set->insert(u);
            }else
            {
                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_map->count(v) && v_edge_set->size() == h - 1 && vertex_map ->at(u) >= k) {
                        --vertex_map->at(u);
                    }
                }
                if (vertex_map->at(u) < k) {
                    evicted_set->insert(u);
                }
            }
        }

        if(evicted_set->size() == vertex_map->size()){
            vertex_map->clear();
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
            for (const auto &u: *evicted_set) {
                vertex_map->erase(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *evicted_set) {
                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_map->count(v)  && v_edge_set->size() >= h && vertex_map->at(v) >= k) {
                        --vertex_map->at(v);
                        if (vertex_map->at(v) < k) {
                            next_evicted_set->insert(v);
                        }
                    }
                }
            }
            removed_set->merge(*evicted_set);
            swap(*evicted_set, *next_evicted_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                                                      uint32_t k,
                                                                                      uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for(const auto &[u, u_degree]:*vertex_map){
            if(vertex_map->at(u) < k){
                evicted_set->insert(u);
            }else if (vertex_edge_size_map->at(u)->count(h - 1)) {
                for (const auto &v: *vertex_edge_size_map->at(u)->at(h - 1)) {
                    if (vertex_map->count(v) && vertex_map->at(u) >= k) {
                        --vertex_map->at(u);
                    }
                }
                if (vertex_map->at(u) < k) {
                    evicted_set->insert(u);
                }
            }
        }

        if(evicted_set->size() == vertex_map->size()){
            vertex_map->clear();
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
            for (const auto &u: *evicted_set) {
                vertex_map->erase(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *evicted_set) {
                auto u_map = vertex_edge_size_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if (vertex_map->count(v) && vertex_map->at(v) >= k) {
                            --vertex_map->at(v);
                            if (vertex_map->at(v) < k) {
                                next_evicted_set->insert(v);
                            }
                        }
                    }
                }
            }
            removed_set->merge(*evicted_set);
            swap(*evicted_set, *next_evicted_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_core(const shared_ptr<temporal_graph>& G,
                                                                                      const shared_ptr<unordered_map<uint32_t,shared_ptr<mutex>>>& vertex_mutex_map,
                                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                      uint32_t k,
                                                                                      uint32_t h,
                                                                                      const shared_ptr<thread_pool>& pool) {
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
                        auto &u = iter->first;
                        if (vertex_degree_map->at(u) < k) {
                            sub_evicted_set->insert(u);
                        } else {
                            for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                                if (vertex_degree_map->count(v) && v_edge_set->size() == h - 1 && vertex_degree_map->at(u) >= k) {
                                    --vertex_degree_map->at(u);
                                }
                            }
                            if (vertex_degree_map->at(u) < k) {
                                sub_evicted_set->insert(u);
                            }
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
            pool->submit_task([=]{
                for (const auto &v: *evicted_set) {
                    vertex_degree_map->erase(v);
                }
            });
            for (const auto &v: *evicted_set) {
                removed_set->insert(v);
            }
            pool->barrier();
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(evicted_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
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

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_core(const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                                                      const shared_ptr<unordered_map<uint32_t,shared_ptr<mutex>>>& vertex_mutex_map,
                                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                      uint32_t k,
                                                                                      uint32_t h,
                                                                                      const shared_ptr<thread_pool>& pool) {
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
                        auto &u = iter->first;
                        if (vertex_degree_map->at(u) < k) {
                            sub_evicted_set->insert(u);
                        } else if (vertex_edge_size_map->at(u)->count(h - 1)) {
                            for (const auto &v: *vertex_edge_size_map->at(u)->at(h - 1)) {
                                if (vertex_degree_map->count(v) && vertex_degree_map->at(u) >= k) {
                                    --vertex_degree_map->at(u);
                                }
                            }
                            if (vertex_degree_map->at(u) < k) {
                                sub_evicted_set->insert(u);
                            }
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
            pool->submit_task([=]{
                for (const auto &v: *evicted_set) {
                    vertex_degree_map->erase(v);
                }
            });
            for (const auto &v: *evicted_set) {
                removed_set->insert(v);
            }
            pool->barrier();
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(evicted_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                        auto &u = *iter1;

                        auto u_map = vertex_edge_size_map->at(u);
                        for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                            for (const auto &v: *iter2->second) {
                                if (vertex_degree_map->count(v) && vertex_degree_map->at(v) >= k) {

                                    vertex_mutex_map->at(v)->lock();
                                    --vertex_degree_map->at(v);
                                    vertex_mutex_map->at(v)->unlock();

                                    if (vertex_degree_map->at(v) < k) {
                                        sub_next_evicted_set->insert(v);
                                    }
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


    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_h_core(const shared_ptr<temporal_graph>& G,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for(const auto &[u, u_degree]:*vertex_map) {
            for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                if (vertex_map->count(v) && v_edge_set->size() == h - 1 && vertex_map->at(u) >= k) {
                    --vertex_map->at(u);
                }
            }
            if (vertex_map->at(u) < k) {
                evicted_set->insert(u);
            }
        }

        if(evicted_set->size() == vertex_map->size()){
            vertex_map->clear();
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
            for (const auto &u: *evicted_set) {
                vertex_map->erase(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *evicted_set) {
                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_map->count(v) && v_edge_set->size() >= h && vertex_map->at(v) >= k) {
                        --vertex_map->at(v);
                        if (vertex_map->at(v) < k) {
                            next_evicted_set->insert(v);
                        }
                    }
                }
            }
            removed_set->merge(*evicted_set);
            swap(*evicted_set, *next_evicted_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_h_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for(const auto &[u, u_degree]:*vertex_map){
            if (vertex_edge_size_map->at(u)->count(h - 1)) {
                for (const auto &v: *vertex_edge_size_map->at(u)->at(h - 1)) {
                    if (vertex_map->count(v) && vertex_map->at(u) >= k) {
                        --vertex_map->at(u);
                    }
                }
                if (vertex_map->at(u) < k) {
                    evicted_set->insert(u);
                }
            }
        }

        if(evicted_set->size() == vertex_map->size()){
            for(const auto&u:*evicted_set){
                vertex_map->erase(u);
            }
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
           for (const auto &u: *evicted_set) {
               vertex_map->erase(u);
           }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *evicted_set) {
                auto u_map = vertex_edge_size_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if (vertex_map->count(v) && vertex_map->at(v) >= k) {
                            --vertex_map->at(v);
                            if (vertex_map->at(v) < k) {
                                next_evicted_set->insert(v);
                            }
                        }
                    }
                }
            }
            removed_set->merge(*evicted_set);
            swap(*evicted_set, *next_evicted_set);
        }
        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_h_core(const shared_ptr<temporal_graph>& G,
                                                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<mutex>>>& vertex_mutex_map,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h,
                                                                                        const shared_ptr<thread_pool>& pool) {
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
                        auto &u = iter->first;
                        for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                            if (vertex_degree_map->count(v) && v_edge_set->size() == h - 1 && vertex_degree_map->at(u) >= k) {
                                --vertex_degree_map->at(u);
                            }
                        }
                        if (vertex_degree_map->at(u) < k) {
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
            pool->submit_task([=]{
                for (const auto &v: *evicted_set) {
                    vertex_degree_map->erase(v);
                }
            });
            for (const auto &v: *evicted_set) {
                removed_set->insert(v);
            }
            pool->barrier();
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(evicted_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
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

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_h_core(const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<mutex>>>& vertex_mutex_map,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h,
                                                                                        const shared_ptr<thread_pool>& pool) {
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
                        auto &u = iter->first;
                        if (vertex_edge_size_map->at(u)->count(h - 1)) {
                            for (const auto &v: *vertex_edge_size_map->at(u)->at(h - 1)) {
                                if (vertex_degree_map->count(v) && vertex_degree_map->at(u) >= k) {
                                    --vertex_degree_map->at(u);
                                }
                            }
                            if (vertex_degree_map->at(u) < k) {
                                sub_evicted_set->insert(u);
                            }
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
            pool->submit_task([=]{
                for (const auto &v: *evicted_set) {
                    vertex_degree_map->erase(v);
                }
            });
            for (const auto &v: *evicted_set) {
                removed_set->insert(v);
            }
            pool->barrier();
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(evicted_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter1 = sub_begin; iter1 != sub_end; ++iter1) {
                        auto &u = *iter1;

                        auto u_map = vertex_edge_size_map->at(u);
                        for (auto iter2 = u_map->lower_bound(h); iter2 != u_map->end(); ++iter2) {
                            for (const auto &v: *iter2->second) {
                                if (vertex_degree_map->count(v) && vertex_degree_map->at(v) >= k) {

                                    vertex_mutex_map->at(v)->lock();
                                    --vertex_degree_map->at(v);
                                    vertex_mutex_map->at(v)->unlock();

                                    if (vertex_degree_map->at(v) < k) {
                                        sub_next_evicted_set->insert(v);
                                    }
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


    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_k_core(const shared_ptr<temporal_graph> &G,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for (const auto&[u,u_degree]: *vertex_map) {
            if (u_degree < k) {
                evicted_set->insert(u);
            }
        }

        if(evicted_set->size() == vertex_map->size()){
            vertex_map->clear();
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
            for (const auto &u: *evicted_set) {
                vertex_map->erase(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &u: *evicted_set) {
                for (const auto &[v, v_edge_set]: *G->get_vertex(u)->get_neighbor_map()) {
                    if (vertex_map->count(v) && v_edge_set->size() >= h && vertex_map->at(v) >= k) {
                        --vertex_map->at(v);
                        if (vertex_map->at(v) < k) {
                            next_evicted_set->insert(v);
                        }
                    }
                }
            }
            removed_set->merge(*evicted_set);
            swap(*evicted_set, *next_evicted_set);
        }
        return removed_set;
    }


    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_k_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h) {
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        for (const auto &[u, u_degree]: *vertex_map) {
            if (vertex_map->at(u) < k) {
                evicted_set->insert(u);
            }
        }

        if(evicted_set->size() == vertex_map->size()){
            vertex_map->clear();
            return evicted_set;
        }

        auto removed_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_set->empty()) {
            for(const auto &u:*evicted_set){
                vertex_map->erase(u);
            }
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            for(const auto &u:*evicted_set) {
                auto u_map = vertex_edge_size_map->at(u);
                for (auto iter = u_map->lower_bound(h); iter != u_map->end(); ++iter) {
                    for (const auto &v: *iter->second) {
                        if (vertex_map->count(v) && vertex_map->at(v) >= k) {
                            --vertex_map->at(v);
                            if (vertex_map->at(v) < k) {
                                next_evicted_set->insert(v);
                            }
                        }
                    }
                }
            }
            removed_set->merge(*evicted_set);
            swap(*evicted_set, *next_evicted_set);
        }

        return removed_set;
    }

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_k_core(const shared_ptr<temporal_graph>& G,
                                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h,
                                                                                        const shared_ptr<thread_pool>& pool) {
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
                        auto &u = iter->first;
                        if (vertex_degree_map->at(u) < k) {
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
            pool->submit_task([=]{
                for(const auto &u:*evicted_set){
                    vertex_degree_map->erase(u);
                }
            });
            for(const auto &u:*evicted_set){
                removed_set->insert(u);
            }
            pool->barrier();
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(evicted_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto sub_begin = *location_vector->at(i);
                    auto sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
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

    shared_ptr<unordered_set<uint32_t>> branch_multiple_core_decomposition::find_k_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t ,shared_ptr<unordered_set<uint32_t>>>>>>& vertex_edge_size_map,
                                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                                        uint32_t k,
                                                                                        uint32_t h,
                                                                                        const shared_ptr<thread_pool>& pool) {
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
                        auto &u = iter->first;
                        if (vertex_degree_map->at(u) < k) {
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
            pool->submit_task([=]{
                for(const auto &u:*evicted_set){
                    vertex_degree_map->erase(u);
                }
            });
            for(const auto &u:*evicted_set){
                removed_set->insert(u);
            };
            pool->barrier();
            auto next_evicted_set = make_shared<unordered_set<uint32_t>>();
            auto location_vector = pool->split_task(evicted_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto sub_next_evicted_set = make_shared<unordered_set<uint32_t>>();

                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);
                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto &u = *iter;

                        auto &u_map = vertex_edge_size_map->at(u);
                        for (auto it = u_map->lower_bound(h); it != u_map->end(); ++it) {
                            for (const auto &v: *it->second) {
                                if (vertex_degree_map->count(v) && vertex_degree_map->at(v) >= k) {

                                    vertex_mutex_map->at(v)->lock();
                                    --vertex_degree_map->at(v);
                                    vertex_mutex_map->at(v)->unlock();

                                    if (vertex_degree_map->at(v) < k) {
                                        sub_next_evicted_set->insert(v);
                                    }
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