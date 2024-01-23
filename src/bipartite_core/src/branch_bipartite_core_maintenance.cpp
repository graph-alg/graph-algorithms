
#include "bipartite_core/branch_bipartite_core_maintenance.h"

namespace scnu{

    void branch_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()) {
            new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
        }
    }

    void branch_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                 const shared_ptr<thread_pool> &pool) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, shared_ptr<bipartite_core_left_store_index>()});
            left_mutex_map->insert({l, shared_ptr<mutex>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()){
            new_right_index_map->insert({r, shared_ptr<bipartite_core_right_store_index>()});
            right_mutex_map->insert({r, shared_ptr<mutex>()});
        }

        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            pool->submit_task([=]{
                new_left_index_map->at(l) = make_shared<bipartite_core_left_store_index>();
                left_mutex_map->at(l) = make_shared<mutex>();
            });
        }
        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            pool->submit_task([=]{
                new_right_index_map->at(r) = make_shared<bipartite_core_right_store_index>();
                right_mutex_map->at(r) = make_shared<mutex>();
            });
        }
        pool->barrier();
    }


    void branch_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                 const shared_ptr<thread_pool>& pool) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, shared_ptr<bipartite_core_left_store_index>()});
            left_core_degree_map->insert({l, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()){
            new_right_index_map->insert({r, shared_ptr<bipartite_core_right_store_index>()});
            right_core_degree_map->insert({r, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }

        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            pool->submit_task([=]{
                new_left_index_map->at(l) = make_shared<bipartite_core_left_store_index>();
                left_core_degree_map->at(l) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                auto degree_map = left_core_degree_map->at(l);
                for(const auto &[r, e]:*l_vertex->get_edge_map()){
                    auto degree = B->get_right_vertex(r)->get_degree();
                    if (!degree_map->count(degree)) {
                        degree_map->insert({degree, make_shared<unordered_set<uint32_t>>()});
                    }
                    degree_map->at(degree)->insert(r);
                }
            });
        }
        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            pool->submit_task([=]{
                new_right_index_map->at(r) = make_shared<bipartite_core_right_store_index>();
                right_core_degree_map->at(r) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                auto degree_map = right_core_degree_map->at(r);
                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    auto degree = B->get_left_vertex(l)->get_degree();
                    if (!degree_map->count(degree)) {
                        degree_map->insert({degree, make_shared<unordered_set<uint32_t>>()});
                    }
                    degree_map->at(degree)->insert(l);
                }
            });
        }
        pool->barrier();
    }

    void branch_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map){

        auto affected_l_set = make_shared<unordered_set<uint32_t>>();
        auto affected_r_set = make_shared<unordered_set<uint32_t>>();

        auto previous_middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto previous_middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();
            B->insert_edge(e);

            affected_l_set->insert(l);
            if (!left_index_map->count(l)) {
                left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});

                previous_middle_inserted_l_map->insert({l, 0});
            }
            left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);

            affected_r_set->insert(r);
            if (!right_index_map->count(r)) {
                right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                previous_middle_inserted_r_map->insert({r, 0});
            }
            right_index_map->at(r)->insert(B->get_right_vertex(r)->get_degree(), 1);
        }


        uint32_t max_k = 1;
        for (uint32_t k = 1; true; ++k) {
            bool flag = false;
            for (const auto &[l, l_index]:*left_index_map) {
                if (l_index->count(k, k)) {
                    flag = true;
                    break;
                }
            }
            if(flag){
                max_k = k;
            }else
            {
                break;
            }
        }

        for(const auto&l :*affected_l_set){
            auto l_vertex = B->get_left_vertex(l);
            auto l_degree = l_vertex->get_degree();
            if (previous_middle_inserted_l_map->count(l)) {
                previous_middle_inserted_l_map->at(l) = l_degree;
            }
            for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                right_index_map->at(r)->insert(1, l_degree);
            }
        }
        for(const auto &r:*affected_r_set){
            auto r_vertex = B->get_right_vertex(r);
            auto r_degree = r_vertex->get_degree();
            if (previous_middle_inserted_r_map->count(r)) {
                previous_middle_inserted_r_map->at(r) = r_degree;
            }
            for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                left_index_map->at(l)->insert(1, r_degree);
            }
        }

        for (uint32_t k = 2; k <= max_k + 1; ++k) {
            {
                auto middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                middle_candidate_graph(B, edge_set, left_index_map, right_index_map,
                                       previous_middle_inserted_l_map, previous_middle_inserted_r_map, k,
                                       middle_inserted_l_map,
                                       middle_inserted_r_map);

                for (const auto &[l, l_degree]: *middle_inserted_l_map) {
                    new_left_index_map->at(l)->insert(k, k);
                }

                for (const auto &[r, r_degree]: *middle_inserted_r_map) {
                    new_right_index_map->at(r)->insert(k, k);
                }

                previous_middle_inserted_l_map->swap(*middle_inserted_l_map);
                previous_middle_inserted_r_map->swap(*middle_inserted_r_map);
            }
            /**
             * @brief update the left path
             */
            {
                auto previous_inserted_l_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_l_map);
                auto previous_inserted_r_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_r_map);

                uint32_t max_i = k;
                for (const auto &[r, r_index]: *right_index_map) {
                    auto value = r_index->get_i(k);

                    if (value > max_i) {
                        max_i = value;
                    }
                }

                for (uint32_t i = k + 1; i <= max_i + 1; ++i) {
                    auto inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                    left_candidate_graph(B, edge_set, left_index_map, right_index_map,
                                         previous_inserted_l_map, previous_inserted_r_map, i, k, inserted_l_map,
                                         inserted_r_map);

                    for (const auto &[l, l_degree]: *inserted_l_map) {
                        new_left_index_map->at(l)->insert(i, k);
                    }

                    for (const auto &[r, r_degree]: *inserted_r_map) {
                        new_right_index_map->at(r)->insert(k, i);
                    }

                    previous_inserted_l_map->swap(*inserted_l_map);
                    previous_inserted_r_map->swap(*inserted_r_map);
                }

                if(!previous_inserted_l_map->empty() && !previous_inserted_r_map->empty()){
                    left_partial_core_decomposition(B, previous_inserted_l_map,
                                                    previous_inserted_r_map,
                                                    new_left_index_map,
                                                    new_right_index_map, k);
                }
            }

            /**
             * @brief update the right path
             */
            {
                auto previous_inserted_l_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_l_map);
                auto previous_inserted_r_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_r_map);

                uint32_t max_j = k;
                for (const auto &[l, l_index]: *left_index_map) {
                    auto value = l_index->get_j(k);
                    if (value > max_j) {
                        max_j = value;
                    }
                }

                for (uint32_t j = k + 1; j <= max_j + 1; ++j) {
                    auto inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                    right_candidate_graph(B, edge_set, left_index_map, right_index_map,
                                          previous_inserted_l_map, previous_inserted_r_map, k, j, inserted_l_map,
                                          inserted_r_map);

                    for (const auto &[l, l_degree]: *inserted_l_map) {
                        new_left_index_map->at(l)->insert(k, j);
                    }

                    for (const auto &[r, r_degree]: *inserted_r_map) {
                        new_right_index_map->at(r)->insert(j, k);
                    }

                    previous_inserted_l_map->swap(*inserted_l_map);
                    previous_inserted_r_map->swap(*inserted_r_map);
                }

                if(!previous_inserted_l_map->empty() && !previous_inserted_r_map->empty()){
                    right_partial_core_decomposition(B, previous_inserted_l_map, previous_inserted_r_map,
                                                     new_left_index_map,
                                                     new_right_index_map, k);
                }
            }
        }

        /**
         * @brief if (max_k + 1)-core is not empty, decompose the remain graph
         */
        if(!previous_middle_inserted_l_map->empty() && !previous_middle_inserted_r_map->empty()){
            middle_partial_core_decomposition(B, previous_middle_inserted_l_map,
                                              previous_middle_inserted_r_map,
                                              new_left_index_map,
                                              new_right_index_map, max_k + 2);
        }


        for(const auto &[l, l_index]:*new_left_index_map){
            for (const auto &[i, j]: *l_index->get_index_map()) {
                left_index_map->at(l)->insert(i, j);
            }
            l_index->clear();
        }

        for(const auto &[r, r_index]:*new_right_index_map){
            for (const auto &[j, i]: *r_index->get_index_map()) {
                right_index_map->at(r)->insert(j, i);
            }
            r_index->clear();
        }
    }

    void branch_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   const shared_ptr<thread_pool> &pool) {
        {
            auto affected_l_set = make_shared<unordered_set<uint32_t>>();
            auto affected_r_set = make_shared<unordered_set<uint32_t>>();

            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                affected_l_set->insert(l);
                affected_r_set->insert(r);

                B->insert_edge(e);

                if (!left_index_map->count(l)) {
                    left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                    new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                }
                left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);


                if (!right_index_map->count(r)) {
                    right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                    new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                }
                right_index_map->at(r)->insert(B->get_right_vertex(r)->get_degree(), 1);
            }

            /**
             * @brief update degree index
             */
            for(const auto&l :*affected_l_set){
                auto l_vertex = B->get_left_vertex(l);
                auto l_degree = l_vertex->get_degree();
                for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                    right_index_map->at(r)->insert(1, l_degree);
                }
            }
            for(const auto &r:*affected_r_set){
                auto r_vertex = B->get_right_vertex(r);
                auto r_degree = r_vertex->get_degree();
                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    left_index_map->at(l)->insert(1, r_degree);
                }
            }
        }

        uint32_t max_k = 1;
        {
            for (uint32_t k = 1; true; ++k) {
                bool flag = false;
                for (const auto &[l, l_index]:*left_index_map) {
                    if (l_index->count(k, k)) {
                        flag = true;
                        break;
                    }
                }
                if(flag){
                    max_k = k;
                }else
                {
                    break;
                }
            }
        }

        auto previous_middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto previous_middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        {
            auto middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            uint32_t k = 2;
            candidate_graph(B, edge_set, left_index_map, right_index_map, k,
                            middle_inserted_l_map,
                            middle_inserted_r_map);

            for (const auto &[l, l_degree]: *middle_inserted_l_map) {
                left_index_map->at(l)->insert(k, k);
            }

            for (const auto &[r, r_degree]: *middle_inserted_r_map) {
                right_index_map->at(r)->insert(k, k);
            }

            previous_middle_inserted_l_map->swap(*middle_inserted_l_map);
            previous_middle_inserted_r_map->swap(*middle_inserted_r_map);
        }

        auto global_left_mutex = make_shared<mutex>();
        auto global_right_mutex = make_shared<mutex>();

        for (uint32_t k = 2; k <= max_k + 1; ++k) {
            /**
             * @brief update the left path
             */
            {
                auto previous_inserted_l_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_l_map);
                auto previous_inserted_r_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_r_map);

                pool->submit_task([=]{
                    uint32_t max_i = k;
                    for (const auto &[r, r_index]: *right_index_map) {
                        auto value = r_index->get_i(k);

                        if (value > max_i) {
                            max_i = value;
                        }
                    }

                    for (uint32_t i = k + 1; i <= max_i + 1; ++i) {
                        auto inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                        left_candidate_graph(B, edge_set, left_index_map, right_index_map,
                                             previous_inserted_l_map, previous_inserted_r_map, i, k, inserted_l_map,
                                             inserted_r_map);

                        global_left_mutex->lock();
                        for (const auto &[l, l_degree]: *inserted_l_map) {
                            new_left_index_map->at(l)->insert(i, k);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for (const auto &[r, r_degree]: *inserted_r_map) {
                            new_right_index_map->at(r)->insert(k, i);
                        }
                        global_right_mutex->unlock();

                        previous_inserted_l_map->swap(*inserted_l_map);
                        previous_inserted_r_map->swap(*inserted_r_map);
                    }

                    if(!previous_inserted_l_map->empty() && !previous_inserted_r_map->empty()){
                        left_partial_core_decomposition(B, global_left_mutex, global_right_mutex, previous_inserted_l_map,
                                                        previous_inserted_r_map,
                                                        new_left_index_map,
                                                        new_right_index_map,  k);
                    }
                });
            }

            /**
             * @brief update the right path
             */
            {
                auto previous_inserted_l_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_l_map);
                auto previous_inserted_r_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_r_map);

                pool->submit_task([=]{
                    uint32_t max_j = k;
                    for(const auto &[l, l_index]:*left_index_map){
                        auto value = l_index->get_j(k);
                        if(value > max_j){
                            max_j = value;
                        }
                    }

                    for (uint32_t j = k + 1; j <= max_j + 1; ++j) {
                        auto inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                        right_candidate_graph(B, edge_set, left_index_map, right_index_map,
                                              previous_inserted_l_map, previous_inserted_r_map, k, j, inserted_l_map,
                                              inserted_r_map);

                        global_left_mutex->lock();
                        for (const auto &[l, l_degree]: *inserted_l_map) {
                            new_left_index_map->at(l)->insert(k, j);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for (const auto &[r, r_degree]: *inserted_r_map) {
                            new_right_index_map->at(r)->insert(j, k);
                        }
                        global_right_mutex->unlock();

                        previous_inserted_l_map->swap(*inserted_l_map);
                        previous_inserted_r_map->swap(*inserted_r_map);
                    }
                    if(!previous_inserted_l_map->empty() && !previous_inserted_r_map->empty()){
                        right_partial_core_decomposition(B, global_left_mutex, global_right_mutex,
                                                         previous_inserted_l_map, previous_inserted_r_map,
                                                         new_left_index_map,
                                                         new_right_index_map, k);
                    }
                });
            }
            {
                auto middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                middle_candidate_graph(B, edge_set, left_index_map, right_index_map,
                                       previous_middle_inserted_l_map, previous_middle_inserted_r_map, k + 1,
                                       middle_inserted_l_map,
                                       middle_inserted_r_map);

                global_left_mutex->lock();
                for (const auto &[l, l_degree]: *middle_inserted_l_map) {
                    new_left_index_map->at(l)->insert(k + 1, k + 1);
                }
                global_left_mutex->unlock();

                global_right_mutex->lock();
                for (const auto &[r, r_degree]: *middle_inserted_r_map) {
                    new_right_index_map->at(r)->insert(k + 1, k + 1);
                }
                global_right_mutex->unlock();

                previous_middle_inserted_l_map->swap(*middle_inserted_l_map);
                previous_middle_inserted_r_map->swap(*middle_inserted_r_map);
            }
        }
        /**
         * @brief if (max_k + 1)-core is not empty, decompose the remain graph
         */
        if (!previous_middle_inserted_l_map->empty() && !previous_middle_inserted_r_map->empty()) {
            middle_partial_core_decomposition(B, global_left_mutex,
                                              global_right_mutex, previous_middle_inserted_l_map,
                                              previous_middle_inserted_r_map,
                                              new_left_index_map,
                                              new_right_index_map, max_k + 1, pool);
        }
        pool->barrier();

        for(const auto &[l, l_index]:*new_left_index_map){
            for(const auto &[i,j]:*l_index->get_index_map()){
                left_index_map->at(l)->insert(i, j);
            }
            l_index->clear();
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            for (const auto &[j, i]: *r_index->get_index_map()) {
                right_index_map->at(r)->insert(j, i);
            }
            r_index->clear();
        }
    }

    void branch_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                   const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                   const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                   const shared_ptr<thread_pool> &pool) {

        auto middle_inserted_l_set = make_shared<unordered_set<uint32_t>>();
        auto middle_inserted_r_set = make_shared<unordered_set<uint32_t>>();
        {
            auto affected_l_set = make_shared<unordered_set<uint32_t>>();
            auto affected_r_set = make_shared<unordered_set<uint32_t>>();

            auto middle_core_rem_degree_map = core_rem_degree_index->get_middle_map();
            auto middle_core_degree_map = core_degree_index->get_middle_map();

            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                affected_l_set->insert(l);
                affected_r_set->insert(r);

                B->insert_edge(e);

                if (!left_index_map->count(l)) {
                    left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                    new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                    middle_core_rem_degree_map->insert({l, 0});
                    middle_core_degree_map->insert({l, 0});

                    middle_inserted_l_set->insert(l);
                    left_mutex_map->insert({l, make_shared<mutex>()});
                }
                left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);

                if (!right_index_map->count(r)) {
                    right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                    new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                    middle_core_rem_degree_map->insert({r, 0});
                    middle_core_degree_map->insert({r, 0});

                    middle_inserted_r_set->insert(r);
                    right_mutex_map->insert({r, make_shared<mutex>()});
                }
                right_index_map->at(r)->insert(B->get_right_vertex(r)->get_degree(), 1);
            }

            /**
             * @brief update degree index
             */

            for (const auto &l: *affected_l_set) {
                auto l_vertex = B->get_left_vertex(l);
                auto l_degree = l_vertex->get_degree();

                if(left_index_map->at(l)->count(1, 1) && !left_index_map->at(l)->count(2, 2)){
                    middle_core_degree_map->at(l) = l_degree;
                }

                for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                    right_index_map->at(r)->insert(1, l_degree);
                }
            }
            for (const auto &r: *affected_r_set) {
                auto r_vertex = B->get_right_vertex(r);
                auto r_degree = r_vertex->get_degree();

                if(right_index_map->at(r)->count(1, 1) && !right_index_map->at(r)->count(2, 2)){
                    middle_core_degree_map->at(r) = r_degree;
                }

                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    left_index_map->at(l)->insert(1, r_degree);
                }
            }
        }

        uint32_t max_k = 1;
        {
            for(const auto &[k, k_list]:*core_order_index->get_middle_map()){
                if (!k_list->empty() && k > max_k) {
                    max_k = k;
                }
            }
        }

        for (uint32_t k = 2; k <= max_k + 1; ++k) {
            auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            {
                auto k_core_order_map = core_order_index->get_middle_map();
                auto k_core_rem_degree_map = core_rem_degree_index->get_middle_map();
                auto k_core_degree_map = core_degree_index->get_middle_map();

                middle_candidate_graph(B, edge_set,
                                       middle_inserted_l_set,
                                       middle_inserted_r_set,
                                       left_index_map,
                                       right_index_map,
                                       k_core_order_map,
                                       k_core_rem_degree_map,
                                       k_core_degree_map,
                                       k);

                for (const auto &l: *middle_inserted_l_set) {
                    left_mutex_map->at(l)->lock();
                    new_left_index_map->at(l)->insert(k, k);
                    left_mutex_map->at(l)->unlock();
                }

                for (const auto &r: *middle_inserted_r_set) {
                    right_mutex_map->at(r)->lock();
                    new_right_index_map->at(r)->insert(k, k);
                    right_mutex_map->at(r)->unlock();
                }

                update_middle_core_degree_map(B, edge_set, left_index_map, right_index_map, middle_inserted_l_set,
                                              middle_inserted_r_set,
                                              k_core_degree_map,
                                              k);

                {
                    for (const auto &l: *middle_inserted_l_set) {
                        vertex_degree_map->insert({l, core_degree_index->get_middle_map()->at(l)});
                    }
                    for (const auto &r: *middle_inserted_r_set) {
                        vertex_degree_map->insert({r, core_degree_index->get_middle_map()->at(r)});
                    }

                    if (k == max_k + 1) {
                        core_order_index->insert_left(k,
                                                      make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                        core_rem_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                                vertex_degree_map->size()));
                        core_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                                vertex_degree_map->size()));

                        core_order_index->insert_right(k,
                                                       make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                        core_rem_degree_index->insert_right(k,
                                                            make_shared<unordered_map<uint32_t, uint32_t>>(
                                                                    vertex_degree_map->size()));
                        core_degree_index->insert_right(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                                vertex_degree_map->size()));
                    }
                }
            }
            pool->submit_task([=] {
                /**
                 * @brief update the left path
                 */
                {
                    auto k_core_order_map = core_order_index->get_left_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_left_map(k);
                    auto k_core_degree_map = core_degree_index->get_left_map(k);

                    auto inserted_l_set = make_shared<unordered_set<uint32_t>>();
                    auto inserted_r_set = make_shared<unordered_set<uint32_t>>();

                    for (const auto &[v, degree]: *vertex_degree_map) {
                        if (B->get_left_vertex(v)) {
                            inserted_l_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        } else {
                            inserted_r_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        }
                    }

                    for (const auto &e: *edge_set) {
                        auto l = e->get_left_vertex_id();
                        auto r = e->get_right_vertex_id();

                        if (left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)
                            && !(left_index_map->at(l)->count(k + 1, k) && right_index_map->at(r)->count(k, k + 1))) {

                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k + 1, k)) {
                                ++k_core_degree_map->at(l);
                            }

                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k, k + 1)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &l: *inserted_l_set) {
                        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k, k + 1)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &r: *inserted_r_set) {
                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k + 1, k)) {
                                ++k_core_degree_map->at(l);
                            }
                        }
                    }

                    uint32_t max_i = k;
                    for (const auto &[i, i_list]: *core_order_index->get_left_map(k)) {
                        if (!i_list->empty() && i > max_i) {
                            max_i = i;
                        }
                    }

                    for (uint32_t i = k + 1; i <= max_i + 1; ++i) {
                        left_candidate_graph(B, edge_set,
                                             inserted_l_set,
                                             inserted_r_set,
                                             left_index_map,
                                             right_index_map,
                                             k_core_order_map,
                                             k_core_rem_degree_map,
                                             k_core_degree_map,
                                             i, k);

                        for (const auto &l: *inserted_l_set) {
                            left_mutex_map->at(l)->lock();
                            new_left_index_map->at(l)->insert(i, k);
                            left_mutex_map->at(l)->unlock();
                        }

                        for (const auto &r: *inserted_r_set) {
                            right_mutex_map->at(r)->lock();
                            new_right_index_map->at(r)->insert(k, i);
                            right_mutex_map->at(r)->unlock();
                        }

                        update_left_core_degree_map(B, edge_set, left_index_map, right_index_map, inserted_l_set,
                                                    inserted_r_set,
                                                    core_degree_index->get_left_map(k), i, k);
                    }

                    if (!inserted_l_set->empty() && !inserted_r_set->empty()) {

                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                inserted_l_set->size() + inserted_r_set->size());
                        for (const auto &l: *inserted_l_set) {
                            sub_vertex_degree_map->insert({l, k_core_degree_map->at(l)});
                        }
                        for (const auto &r: *inserted_r_set) {
                            sub_vertex_degree_map->insert({r, k_core_degree_map->at(r)});
                        }


                        left_partial_core_decomposition(B,
                                                        left_mutex_map,
                                                        right_mutex_map,
                                                        inserted_l_set,
                                                        inserted_r_set,
                                                        sub_vertex_degree_map,
                                                        new_left_index_map,
                                                        new_right_index_map,
                                                        k_core_order_map,
                                                        k_core_rem_degree_map,
                                                        k_core_degree_map,
                                                        k);
                    }
                }

                /**
                 * @brief right path
                 */
                {
                    auto k_core_order_map = core_order_index->get_right_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_right_map(k);
                    auto k_core_degree_map = core_degree_index->get_right_map(k);

                    auto inserted_l_set = make_shared<unordered_set<uint32_t>>();
                    auto inserted_r_set = make_shared<unordered_set<uint32_t>>();

                    for (const auto &[v, degree]: *vertex_degree_map) {
                        if (B->get_left_vertex(v)) {
                            inserted_l_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        } else {
                            inserted_r_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        }
                    }

                    for (const auto &e: *edge_set) {
                        auto l = e->get_left_vertex_id();
                        auto r = e->get_right_vertex_id();

                        if (left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)
                            && !(left_index_map->at(l)->count(k, k + 1) && right_index_map->at(r)->count(k + 1, k))) {

                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k, k + 1)) {
                                ++k_core_degree_map->at(l);
                            }

                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k + 1, k)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &l: *inserted_l_set) {
                        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k + 1, k)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &r: *inserted_r_set) {
                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k, k + 1)) {
                                ++k_core_degree_map->at(l);
                            }
                        }
                    }

                    uint32_t max_j = k;
                    for (const auto &[j, j_list]: *core_order_index->get_right_map(k)) {
                        if (!j_list->empty() && j > max_j) {
                            max_j = j;
                        }
                    }

                    for (uint32_t j = k + 1; j <= max_j + 1; ++j) {
                        right_candidate_graph(B, edge_set,
                                              inserted_l_set,
                                              inserted_r_set,
                                              left_index_map,
                                              right_index_map,
                                              k_core_order_map,
                                              k_core_rem_degree_map,
                                              k_core_degree_map,
                                              k, j);

                        for (const auto &l: *inserted_l_set) {
                            left_mutex_map->at(l)->lock();
                            new_left_index_map->at(l)->insert(k, j);
                            left_mutex_map->at(l)->unlock();
                        }

                        for (const auto &r: *inserted_r_set) {
                            right_mutex_map->at(r)->lock();
                            new_right_index_map->at(r)->insert(j, k);
                            right_mutex_map->at(r)->unlock();
                        }

                        update_right_core_degree_map(B, edge_set, left_index_map, right_index_map, inserted_l_set,
                                                     inserted_r_set,
                                                     k_core_degree_map, k, j);

                    }
                    if (!inserted_l_set->empty() && !inserted_r_set->empty()) {
                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                inserted_l_set->size() + inserted_r_set->size());
                        for (const auto &l: *inserted_l_set) {
                            sub_vertex_degree_map->insert({l, k_core_degree_map->at(l)});
                        }
                        for (const auto &r: *inserted_r_set) {
                            sub_vertex_degree_map->insert({r, k_core_degree_map->at(r)});
                        }
                        right_partial_core_decomposition(B,
                                                         left_mutex_map,
                                                         right_mutex_map,
                                                         inserted_l_set,
                                                         inserted_r_set,
                                                         sub_vertex_degree_map,
                                                         new_left_index_map,
                                                         new_right_index_map,
                                                         k_core_order_map,
                                                         k_core_rem_degree_map,
                                                         k_core_degree_map,
                                                         k);
                    }
                }
            });
        }

        /**
         * @brief if (max_k + 1)-core is not empty, decompose the remain graph
         */
        if (!middle_inserted_l_set->empty() && !middle_inserted_r_set->empty()) {
            middle_partial_core_decomposition(B, left_mutex_map,
                                              right_mutex_map,
                                              middle_inserted_l_set,
                                              middle_inserted_r_set,
                                              new_left_index_map,
                                              new_right_index_map,
                                              core_order_index,
                                              core_rem_degree_index,
                                              core_degree_index,
                                              max_k + 1,
                                              pool);
        }

        pool->barrier();

        for (const auto &[l, l_index]: *new_left_index_map) {
            pool->submit_task([=]{
                for (const auto &[i, j]: *l_index->get_index_map()) {
                    left_index_map->at(l)->insert(i, j);
                }
                l_index->clear();
            });
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            pool->submit_task([=] {
                for (const auto &[j, i]: *r_index->get_index_map()) {
                    right_index_map->at(r)->insert(j, i);
                }
                r_index->clear();
            });
        }

        pool->barrier();
    }

    void branch_bipartite_core_maintenance::insert2(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                    const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                    const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                    const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                    const shared_ptr<thread_pool> &pool) {

        auto middle_inserted_l_set = make_shared<unordered_set<uint32_t>>();
        auto middle_inserted_r_set = make_shared<unordered_set<uint32_t>>();
        {
            auto affected_l_set = make_shared<unordered_set<uint32_t>>();
            auto affected_r_set = make_shared<unordered_set<uint32_t>>();

            auto middle_core_rem_degree_map = core_rem_degree_index->get_middle_map();
            auto middle_core_degree_map = core_degree_index->get_middle_map();

            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                affected_l_set->insert(l);
                affected_r_set->insert(r);

                B->insert_edge(e);

                if (!left_index_map->count(l)) {
                    left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                    new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                    middle_core_rem_degree_map->insert({l, 0});
                    middle_core_degree_map->insert({l, 0});

                    middle_inserted_l_set->insert(l);
                    left_mutex_map->insert({l, make_shared<mutex>()});
                }
                left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);

                if (!right_index_map->count(r)) {
                    right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                    new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                    middle_core_rem_degree_map->insert({r, 0});
                    middle_core_degree_map->insert({r, 0});

                    middle_inserted_r_set->insert(r);
                    right_mutex_map->insert({r, make_shared<mutex>()});
                }
                right_index_map->at(r)->insert(B->get_right_vertex(r)->get_degree(), 1);
            }

            /**
             * @brief update degree index
             */

            for (const auto &l: *affected_l_set) {
                auto l_vertex = B->get_left_vertex(l);
                auto l_degree = l_vertex->get_degree();

                if (left_index_map->at(l)->count(1, 1) && !left_index_map->at(l)->count(2, 2)) {
                    middle_core_degree_map->at(l) = l_degree;
                }

                for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                    right_index_map->at(r)->insert(1, l_degree);
                }
            }
            for (const auto &r: *affected_r_set) {
                auto r_vertex = B->get_right_vertex(r);
                auto r_degree = r_vertex->get_degree();

                if (right_index_map->at(r)->count(1, 1) && !right_index_map->at(r)->count(2, 2)) {
                    middle_core_degree_map->at(r) = r_degree;
                }

                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    left_index_map->at(l)->insert(1, r_degree);
                }
            }
        }

        uint32_t max_k = 1;
        {
            for (const auto &[k, k_list]: *core_order_index->get_middle_map()) {
                if (!k_list->empty() && k > max_k) {
                    max_k = k;
                }
            }
        }


        auto k_vector = make_shared<vector<shared_ptr<unordered_map<uint32_t, uint32_t>>>>(max_k + 2);
        for (uint32_t k = 2; k <= max_k + 1; ++k) {
            auto k_core_order_map = core_order_index->get_middle_map();
            auto k_core_rem_degree_map = core_rem_degree_index->get_middle_map();
            auto k_core_degree_map = core_degree_index->get_middle_map();

            middle_candidate_graph(B, edge_set,
                                   middle_inserted_l_set,
                                   middle_inserted_r_set,
                                   left_index_map,
                                   right_index_map,
                                   k_core_order_map,
                                   k_core_rem_degree_map,
                                   k_core_degree_map,
                                   k);

            for (const auto &l: *middle_inserted_l_set) {
                new_left_index_map->at(l)->insert(k, k);
            }

            for (const auto &r: *middle_inserted_r_set) {
                new_right_index_map->at(r)->insert(k, k);
            }

            update_middle_core_degree_map(B, edge_set, left_index_map, right_index_map, middle_inserted_l_set,
                                          middle_inserted_r_set,
                                          k_core_degree_map,
                                          k);

            {
                k_vector->at(k) = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto vertex_degree_map = k_vector->at(k);
                for (const auto &l: *middle_inserted_l_set) {
                    vertex_degree_map->insert({l, core_degree_index->get_middle_map()->at(l)});
                }
                for (const auto &r: *middle_inserted_r_set) {
                    vertex_degree_map->insert({r, core_degree_index->get_middle_map()->at(r)});
                }

                if (k == max_k + 1) {
                    core_order_index->insert_left(k,
                                                  make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                    core_rem_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                            vertex_degree_map->size()));
                    core_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                            vertex_degree_map->size()));

                    core_order_index->insert_right(k,
                                                   make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                    core_rem_degree_index->insert_right(k,
                                                        make_shared<unordered_map<uint32_t, uint32_t>>(
                                                                vertex_degree_map->size()));
                    core_degree_index->insert_right(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                            vertex_degree_map->size()));
                }
            }
        }


        for (uint32_t k = 2; k <= max_k + 1; ++k) {
            pool->submit_task([=] {
                auto vertex_degree_map = k_vector->at(k);
                /**
                 * @brief update the left path
                 */
                {
                    auto k_core_order_map = core_order_index->get_left_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_left_map(k);
                    auto k_core_degree_map = core_degree_index->get_left_map(k);

                    auto inserted_l_set = make_shared<unordered_set<uint32_t>>();
                    auto inserted_r_set = make_shared<unordered_set<uint32_t>>();

                    for (const auto &[v, degree]: *vertex_degree_map) {
                        if (B->get_left_vertex(v)) {
                            inserted_l_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        } else {
                            inserted_r_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        }
                    }

                    for (const auto &e: *edge_set) {
                        auto l = e->get_left_vertex_id();
                        auto r = e->get_right_vertex_id();

                        if (left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)
                            && !(left_index_map->at(l)->count(k + 1, k) && right_index_map->at(r)->count(k, k + 1))) {

                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k + 1, k)) {
                                ++k_core_degree_map->at(l);
                            }

                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k, k + 1)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &l: *inserted_l_set) {
                        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k, k + 1)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &r: *inserted_r_set) {
                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k + 1, k)) {
                                ++k_core_degree_map->at(l);
                            }
                        }
                    }

                    uint32_t max_i = k;
                    for (const auto &[i, i_list]: *core_order_index->get_left_map(k)) {
                        if (!i_list->empty() && i > max_i) {
                            max_i = i;
                        }
                    }

                    for (uint32_t i = k + 1; i <= max_i + 1; ++i) {
                        left_candidate_graph(B, edge_set,
                                             inserted_l_set,
                                             inserted_r_set,
                                             left_index_map,
                                             right_index_map,
                                             k_core_order_map,
                                             k_core_rem_degree_map,
                                             k_core_degree_map,
                                             i, k);

                        for (const auto &l: *inserted_l_set) {
                            left_mutex_map->at(l)->lock();
                            new_left_index_map->at(l)->insert(i, k);
                            left_mutex_map->at(l)->unlock();
                        }

                        for (const auto &r: *inserted_r_set) {
                            right_mutex_map->at(r)->lock();
                            new_right_index_map->at(r)->insert(k, i);
                            right_mutex_map->at(r)->unlock();
                        }

                        update_left_core_degree_map(B, edge_set, left_index_map, right_index_map, inserted_l_set,
                                                    inserted_r_set,
                                                    core_degree_index->get_left_map(k), i, k);
                    }

                    if (!inserted_l_set->empty() && !inserted_r_set->empty()) {

                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                inserted_l_set->size() + inserted_r_set->size());
                        for (const auto &l: *inserted_l_set) {
                            sub_vertex_degree_map->insert({l, k_core_degree_map->at(l)});
                        }
                        for (const auto &r: *inserted_r_set) {
                            sub_vertex_degree_map->insert({r, k_core_degree_map->at(r)});
                        }


                        left_partial_core_decomposition(B,
                                                        left_mutex_map,
                                                        right_mutex_map,
                                                        inserted_l_set,
                                                        inserted_r_set,
                                                        sub_vertex_degree_map,
                                                        new_left_index_map,
                                                        new_right_index_map,
                                                        k_core_order_map,
                                                        k_core_rem_degree_map,
                                                        k_core_degree_map,
                                                        k);
                    }
                }
            });
        }

        for (uint32_t k = 2; k <= max_k + 1; ++k) {
            pool->submit_task([=] {
                auto vertex_degree_map = k_vector->at(k);
                /**
                 * @brief right path
                 */
                {
                    auto k_core_order_map = core_order_index->get_right_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_right_map(k);
                    auto k_core_degree_map = core_degree_index->get_right_map(k);

                    auto inserted_l_set = make_shared<unordered_set<uint32_t>>();
                    auto inserted_r_set = make_shared<unordered_set<uint32_t>>();

                    for (const auto &[v, degree]: *vertex_degree_map) {
                        if (B->get_left_vertex(v)) {
                            inserted_l_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        } else {
                            inserted_r_set->insert(v);
                            k_core_degree_map->insert({v, degree});
                            k_core_rem_degree_map->insert({v, 0});
                        }
                    }

                    for (const auto &e: *edge_set) {
                        auto l = e->get_left_vertex_id();
                        auto r = e->get_right_vertex_id();

                        if (left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)
                            && !(left_index_map->at(l)->count(k, k + 1) && right_index_map->at(r)->count(k + 1, k))) {

                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k, k + 1)) {
                                ++k_core_degree_map->at(l);
                            }

                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k + 1, k)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &l: *inserted_l_set) {
                        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k + 1, k)) {
                                ++k_core_degree_map->at(r);
                            }
                        }
                    }

                    for (const auto &r: *inserted_r_set) {
                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k, k + 1)) {
                                ++k_core_degree_map->at(l);
                            }
                        }
                    }

                    uint32_t max_j = k;
                    for (const auto &[j, j_list]: *core_order_index->get_right_map(k)) {
                        if (!j_list->empty() && j > max_j) {
                            max_j = j;
                        }
                    }

                    for (uint32_t j = k + 1; j <= max_j + 1; ++j) {
                        right_candidate_graph(B, edge_set,
                                              inserted_l_set,
                                              inserted_r_set,
                                              left_index_map,
                                              right_index_map,
                                              k_core_order_map,
                                              k_core_rem_degree_map,
                                              k_core_degree_map,
                                              k, j);

                        for (const auto &l: *inserted_l_set) {
                            left_mutex_map->at(l)->lock();
                            new_left_index_map->at(l)->insert(k, j);
                            left_mutex_map->at(l)->unlock();
                        }

                        for (const auto &r: *inserted_r_set) {
                            right_mutex_map->at(r)->lock();
                            new_right_index_map->at(r)->insert(j, k);
                            right_mutex_map->at(r)->unlock();
                        }

                        update_right_core_degree_map(B, edge_set, left_index_map, right_index_map, inserted_l_set,
                                                     inserted_r_set,
                                                     k_core_degree_map, k, j);

                    }
                    if (!inserted_l_set->empty() && !inserted_r_set->empty()) {
                        auto sub_vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                                inserted_l_set->size() + inserted_r_set->size());
                        for (const auto &l: *inserted_l_set) {
                            sub_vertex_degree_map->insert({l, k_core_degree_map->at(l)});
                        }
                        for (const auto &r: *inserted_r_set) {
                            sub_vertex_degree_map->insert({r, k_core_degree_map->at(r)});
                        }
                        right_partial_core_decomposition(B,
                                                         left_mutex_map,
                                                         right_mutex_map,
                                                         inserted_l_set,
                                                         inserted_r_set,
                                                         sub_vertex_degree_map,
                                                         new_left_index_map,
                                                         new_right_index_map,
                                                         k_core_order_map,
                                                         k_core_rem_degree_map,
                                                         k_core_degree_map,
                                                         k);
                    }
                }
            });
        }

        /**
         * @brief if (max_k + 1)-core is not empty, decompose the remain graph
         */
        if (!middle_inserted_l_set->empty() && !middle_inserted_r_set->empty()) {
            middle_partial_core_decomposition2(B, left_mutex_map,
                                               right_mutex_map,
                                               middle_inserted_l_set,
                                               middle_inserted_r_set,
                                               new_left_index_map,
                                               new_right_index_map,
                                               core_order_index,
                                               core_rem_degree_index,
                                               core_degree_index,
                                               max_k + 1,
                                               pool);
        }

        pool->barrier();

        for (const auto &[l, l_index]: *new_left_index_map) {
            pool->submit_task([=] {
                for (const auto &[i, j]: *l_index->get_index_map()) {
                    left_index_map->at(l)->insert(i, j);
                }
                l_index->clear();
            });
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            pool->submit_task([=] {
                for (const auto &[j, i]: *r_index->get_index_map()) {
                    right_index_map->at(r)->insert(j, i);
                }
                r_index->clear();
            });
        }

        pool->barrier();
    }


    void branch_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   const shared_ptr<thread_pool>& pool){
        {
            auto affected_l_set = make_shared<unordered_set<uint32_t>>();
            auto affected_r_set = make_shared<unordered_set<uint32_t>>();

            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                affected_l_set->insert(l);
                affected_r_set->insert(r);

                B->insert_edge(e);

                if (!left_index_map->count(l)) {
                    left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                    new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                }
                left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);


                if (!right_index_map->count(r)) {
                    right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                    new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                }
                right_index_map->at(r)->insert(B->get_right_vertex(r)->get_degree(), 1);
            }

            /**
             * @brief update degree index
             */
            for(const auto&l :*affected_l_set){
                auto l_vertex = B->get_left_vertex(l);
                auto l_degree = l_vertex->get_degree();

                for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                    right_index_map->at(r)->insert(1, l_degree);
                }
            }
            for(const auto &r:*affected_r_set){
                auto r_vertex = B->get_right_vertex(r);
                auto r_degree = r_vertex->get_degree();

                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    left_index_map->at(l)->insert(1, r_degree);
                }
            }

            update_index_for_insertion(B, edge_set, left_core_degree_map, right_core_degree_map, affected_l_set, affected_r_set);
        }

        uint32_t max_k = 1;
        {
            for (uint32_t k = 1; true; ++k) {
                bool flag = false;
                for (const auto &[l, l_index]:*left_index_map) {
                    if (l_index->count(k, k)) {
                        flag = true;
                        break;
                    }
                }
                if(flag){
                    max_k = k;
                }else
                {
                    break;
                }
            }
        }

        auto global_left_mutex = make_shared<mutex>();
        auto global_right_mutex = make_shared<mutex>();
        auto previous_middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto previous_middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        {
            auto middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            auto middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

            uint32_t k = 2;
            candidate_graph(B, edge_set, left_index_map, right_index_map, k,
                            middle_inserted_l_map,
                            middle_inserted_r_map);


            for (const auto &[l, l_degree]: *middle_inserted_l_map) {
                left_index_map->at(l)->insert(k, k);
            }

            for (const auto &[r, r_degree]: *middle_inserted_r_map) {
                right_index_map->at(r)->insert(k, k);
            }

            previous_middle_inserted_l_map->swap(*middle_inserted_l_map);
            previous_middle_inserted_r_map->swap(*middle_inserted_r_map);
        }



        for (uint32_t k = 2; k <= max_k + 1; ++k) {
            /**
             * @brief update the left path
             */
            {
                auto previous_inserted_l_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_l_map);
                auto previous_inserted_r_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_r_map);

                pool->submit_task([=]{
                    uint32_t max_i = k;
                    for (const auto &[r, r_index]: *right_index_map) {
                        auto value = r_index->get_i(k);

                        if (value > max_i) {
                            max_i = value;
                        }
                    }

                    for (uint32_t i = k + 1; i <= max_i + 1; ++i) {
                        auto inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                        left_candidate_graph(left_core_degree_map, right_core_degree_map, edge_set, left_index_map, right_index_map,
                                             previous_inserted_l_map, previous_inserted_r_map, i, k, inserted_l_map,
                                             inserted_r_map);

                        global_left_mutex->lock();
                        for (const auto &[l, l_degree]: *inserted_l_map) {
                            new_left_index_map->at(l)->insert(i, k);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for (const auto &[r, r_degree]: *inserted_r_map) {
                            new_right_index_map->at(r)->insert(k, i);
                        }
                        global_right_mutex->unlock();

                        previous_inserted_l_map->swap(*inserted_l_map);
                        previous_inserted_r_map->swap(*inserted_r_map);
                    }

                    if(!previous_inserted_l_map->empty() && !previous_inserted_r_map->empty()){
                        left_partial_core_decomposition(left_core_degree_map, right_core_degree_map, global_left_mutex, global_right_mutex, previous_inserted_l_map,
                                                        previous_inserted_r_map,
                                                        new_left_index_map,
                                                        new_right_index_map, k);
                    }
                });
            }

            /**
             * @brief update the right path
             */
            {
                auto previous_inserted_l_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_l_map);
                auto previous_inserted_r_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                        previous_middle_inserted_r_map);

                pool->submit_task([=]{
                    uint32_t max_j = k;
                    for(const auto &[l, l_index]:*left_index_map){
                        auto value = l_index->get_j(k);
                        if(value > max_j){
                            max_j = value;
                        }
                    }

                    for (uint32_t j = k + 1; j <= max_j + 1; ++j) {
                        auto inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                        auto inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                        right_candidate_graph(left_core_degree_map, right_core_degree_map, edge_set, left_index_map, right_index_map,
                                              previous_inserted_l_map, previous_inserted_r_map, k, j, inserted_l_map,
                                              inserted_r_map);

                        global_left_mutex->lock();
                        for (const auto &[l, l_degree]: *inserted_l_map) {
                            new_left_index_map->at(l)->insert(k, j);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for (const auto &[r, r_degree]: *inserted_r_map) {
                            new_right_index_map->at(r)->insert(j, k);
                        }
                        global_right_mutex->unlock();

                        previous_inserted_l_map->swap(*inserted_l_map);
                        previous_inserted_r_map->swap(*inserted_r_map);
                    }
                    if(!previous_inserted_l_map->empty() && !previous_inserted_r_map->empty()){
                        right_partial_core_decomposition(left_core_degree_map, right_core_degree_map, global_left_mutex, global_right_mutex,
                                                         previous_inserted_l_map, previous_inserted_r_map,
                                                         new_left_index_map,
                                                         new_right_index_map, k);
                    }
                });
            }

            {
                auto middle_inserted_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto middle_inserted_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                middle_candidate_graph(left_core_degree_map, right_core_degree_map, edge_set, left_index_map, right_index_map,
                                       previous_middle_inserted_l_map,
                                       previous_middle_inserted_r_map, k + 1,
                                       middle_inserted_l_map,
                                       middle_inserted_r_map);


                global_left_mutex->lock();
                for (const auto &[l, l_degree]: *middle_inserted_l_map) {
                    new_left_index_map->at(l)->insert(k + 1, k + 1);
                }
                global_left_mutex->unlock();

                global_right_mutex->lock();
                for (const auto &[r, r_degree]: *middle_inserted_r_map) {
                    new_right_index_map->at(r)->insert(k + 1, k + 1);
                }
                global_right_mutex->unlock();

                previous_middle_inserted_l_map->clear();
                previous_middle_inserted_r_map->clear();

                previous_middle_inserted_l_map->swap(*middle_inserted_l_map);
                previous_middle_inserted_r_map->swap(*middle_inserted_r_map);
            }
        }
        /**
         * @brief if (max_k + 1)-core is not empty, decompose the remain graph
         */
        if (!previous_middle_inserted_l_map->empty() && !previous_middle_inserted_r_map->empty()) {
            middle_partial_core_decomposition(left_core_degree_map, right_core_degree_map, global_left_mutex,
                                              global_right_mutex, previous_middle_inserted_l_map,
                                              previous_middle_inserted_r_map,
                                              new_left_index_map,
                                              new_right_index_map, max_k + 1, pool);
        }
        pool->barrier();

        for(const auto &[l, l_index]:*new_left_index_map){
            for(const auto &[i,j]:*l_index->get_index_map()){
                left_index_map->at(l)->insert(i, j);
            }
            l_index->clear();
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            for (const auto &[j, i]: *r_index->get_index_map()) {
                right_index_map->at(r)->insert(j, i);
            }
            r_index->clear();
        }
    }



    void branch_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map) {
        auto isolated_left_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_right_vertex_set = make_shared<unordered_set<uint32_t>>();

        B->remove_edge_collection(edge_set, isolated_left_vertex_set, isolated_right_vertex_set);

        uint32_t max_k = 0;
        while(true){
            auto flag = false;
            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                if(left_index_map->at(l)->count(max_k + 1, max_k + 1) && right_index_map->at(r)->count(max_k + 1, max_k + 1)){
                    flag = true;
                    break;
                }
            }
            if(flag){
                ++max_k;
            }else{
                break;
            }
        }

        auto previous_middle_removed_l_set =  make_shared<unordered_set<uint32_t>>();
        auto previous_middle_removed_r_set =  make_shared<unordered_set<uint32_t>>();
        for(uint32_t k = max_k; k >= 2; --k){
            {
                uint32_t max_i = k;
                while (true) {
                    auto flag = false;
                    for (const auto &e: *edge_set) {
                        auto l = e->get_left_vertex_id();
                        auto r = e->get_right_vertex_id();

                        if (left_index_map->at(l)->count(max_i + 1, k) &&
                            right_index_map->at(r)->count(k, max_i + 1)) {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        ++max_i;
                    } else {
                        break;
                    }
                }

                auto previous_removed_l_set =  make_shared<unordered_set<uint32_t>>();
                auto previous_removed_r_set =  make_shared<unordered_set<uint32_t>>();

                for (uint32_t  i = max_i; i > k; --i){
                    auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                    auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                    left_removal_partial_core(B, edge_set, left_index_map, right_index_map, previous_removed_l_set, previous_removed_r_set,i, k, removed_l_set, removed_r_set);

                    for(const auto&l:*removed_l_set){
                        new_left_index_map->at(l)->remove(i, k);
                    }

                    for(const auto&r:*removed_r_set){
                        new_right_index_map->at(r)->remove(k, i);
                    }

                    swap(*previous_removed_l_set, *removed_l_set);
                    swap(*previous_removed_r_set, *removed_r_set);
                }
            }

            {
                uint32_t max_j = k;
                while(true){
                    auto flag = false;
                    for (const auto &e: *edge_set) {
                        auto l = e->get_left_vertex_id();
                        auto r = e->get_right_vertex_id();

                        if (left_index_map->at(l)->count(k, max_j + 1) &&
                            right_index_map->at(r)->count(max_j + 1, k)) {
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        ++max_j;
                    } else {
                        break;
                    }
                }

                auto previous_removed_l_set =  make_shared<unordered_set<uint32_t>>();
                auto previous_removed_r_set =  make_shared<unordered_set<uint32_t>>();

                for (uint32_t  j = max_j; j > k; --j){
                    auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                    auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                    right_removal_partial_core(B, edge_set, left_index_map, right_index_map, previous_removed_l_set, previous_removed_r_set,k, j, removed_l_set, removed_r_set);

                    for(const auto &l:*removed_l_set){
                        new_left_index_map->at(l)->remove(k, j);
                    }

                    for(const auto &r:*removed_r_set){
                        new_right_index_map->at(r)->remove(j, k);
                    }

                    swap(*previous_removed_l_set, *removed_l_set);
                    swap(*previous_removed_r_set, *removed_r_set);
                }
            }

            {

                auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                middle_removal_partial_core(B, edge_set, left_index_map, right_index_map,
                                            previous_middle_removed_l_set, previous_middle_removed_r_set, k, removed_l_set, removed_r_set);

                for(const auto &l:*removed_l_set){
                    new_left_index_map->at(l)->remove(k, k);
                }

                for(const auto &r:*removed_r_set){
                    new_right_index_map->at(r)->remove(k, k);
                }

                swap(*previous_middle_removed_l_set, *removed_l_set);
                swap(*previous_middle_removed_r_set, *removed_r_set);
            }
        }

        for(const auto&[l, l_index]:*new_left_index_map) {
            for (const auto &[i, j]: *l_index->get_index_map()) {
                left_index_map->at(l)->remove(i, j - 1);
                if(left_index_map->at(l)->get_j(i) == 0){
                    left_index_map->at(l)->remove(i);
                }
            }
            l_index->clear();
        }

        for(const auto &[r, r_index]:*new_right_index_map){
            for(const auto [j, i]:*r_index->get_index_map()){
                right_index_map->at(r)->remove(j, i - 1);
                if( right_index_map->at(r)->get_i(j) == 0){
                    right_index_map->at(r)->remove(j);
                }
            }
            r_index->clear();
        }

        update_trivial_bipartite_cores(B, edge_set, left_index_map, right_index_map);

        for(const auto&l:*isolated_left_vertex_set){
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }

        for(const auto &r:*isolated_right_vertex_set){
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }
    }

    void branch_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   const shared_ptr<thread_pool>& pool) {
        uint32_t max_k = 0;
        while(true){
            auto flag = false;
            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                if(left_index_map->at(l)->count(max_k + 1, max_k + 1) && right_index_map->at(r)->count(max_k + 1, max_k + 1)){
                    flag = true;
                    break;
                }
            }
            if(flag){
                ++max_k;
            }else{
                break;
            }
        }

        auto isolated_left_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_right_vertex_set = make_shared<unordered_set<uint32_t>>();

        B->remove_edge_collection(edge_set, isolated_left_vertex_set, isolated_right_vertex_set);

        auto global_left_mutex = make_shared<mutex>();
        auto global_right_mutex = make_shared<mutex>();

        for(uint32_t k = max_k; k >= 2; --k){
            {
                pool->submit_task([=]{
                    uint32_t max_i = k;
                    while (true) {
                        auto flag = false;
                        for (const auto &e: *edge_set) {
                            auto l = e->get_left_vertex_id();
                            auto r = e->get_right_vertex_id();

                            if (left_index_map->at(l)->count(max_i + 1, k) &&
                                right_index_map->at(r)->count(k, max_i + 1)) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) {
                            ++max_i;
                        } else {
                            break;
                        }
                    }

                    auto previous_removed_l_set =  make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set =  make_shared<unordered_set<uint32_t>>();

                    for (uint32_t  i = max_i; i > k; --i){
                        auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                        auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                        left_removal_partial_core(B, edge_set, left_index_map, right_index_map, previous_removed_l_set, previous_removed_r_set,i, k, removed_l_set, removed_r_set);

                        global_left_mutex->lock();
                        for(const auto&l:*removed_l_set){
                            new_left_index_map->at(l)->remove(i, k);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for(const auto&r:*removed_r_set){
                            new_right_index_map->at(r)->remove(k, i);
                        }
                        global_right_mutex->unlock();

                        swap(*previous_removed_l_set, *removed_l_set);
                        swap(*previous_removed_r_set, *removed_r_set);
                    }
                });
            }

            {
                pool->submit_task([=]{
                    uint32_t max_j = k;
                    while(true){
                        auto flag = false;
                        for (const auto &e: *edge_set) {
                            auto l = e->get_left_vertex_id();
                            auto r = e->get_right_vertex_id();

                            if (left_index_map->at(l)->count(k, max_j + 1) &&
                                right_index_map->at(r)->count(max_j + 1, k)) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) {
                            ++max_j;
                        } else {
                            break;
                        }
                    }

                    auto previous_removed_l_set =  make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set =  make_shared<unordered_set<uint32_t>>();

                    for (uint32_t  j = max_j; j > k; --j){
                        auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                        auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                        right_removal_partial_core(B, edge_set, left_index_map, right_index_map, previous_removed_l_set, previous_removed_r_set,k, j, removed_l_set, removed_r_set);

                        global_left_mutex->lock();
                        for(const auto &l:*removed_l_set){
                            new_left_index_map->at(l)->remove(k, j);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for(const auto &r:*removed_r_set){
                            new_right_index_map->at(r)->remove(j, k);
                        }
                        global_right_mutex->unlock();

                        swap(*previous_removed_l_set, *removed_l_set);
                        swap(*previous_removed_r_set, *removed_r_set);
                    }
                });
            }
        }
        auto previous_middle_removed_l_set =  make_shared<unordered_set<uint32_t>>();
        auto previous_middle_removed_r_set =  make_shared<unordered_set<uint32_t>>();

        for(uint32_t k = max_k; k >= 2; --k){
            auto removed_l_set = make_shared<unordered_set<uint32_t>>();
            auto removed_r_set = make_shared<unordered_set<uint32_t>>();

            middle_removal_partial_core(B, edge_set, left_index_map, right_index_map,
                                        previous_middle_removed_l_set, previous_middle_removed_r_set, k, removed_l_set, removed_r_set);

            global_left_mutex->lock();
            for(const auto &l:*removed_l_set){
                new_left_index_map->at(l)->remove(k, k);
            }
            global_left_mutex->unlock();

            global_right_mutex->lock();
            for(const auto &r:*removed_r_set){
                new_right_index_map->at(r)->remove(k, k);
            }
            global_right_mutex->unlock();

            swap(*previous_middle_removed_l_set, *removed_l_set);
            swap(*previous_middle_removed_r_set, *removed_r_set);
        }

        pool->barrier();

        for(const auto&[l, l_index]:*new_left_index_map) {
            for (const auto &[i, j]: *l_index->get_index_map()) {
                left_index_map->at(l)->remove(i, j - 1);
                if(left_index_map->at(l)->get_j(i) == 0){
                    left_index_map->at(l)->remove(i);
                }
            }
            l_index->clear();
        }

        for(const auto &[r, r_index]:*new_right_index_map){
            for(const auto [j, i]:*r_index->get_index_map()){
                right_index_map->at(r)->remove(j, i - 1);
                if( right_index_map->at(r)->get_i(j) == 0){
                    right_index_map->at(r)->remove(j);
                }
            }
            r_index->clear();
        }

        update_trivial_bipartite_cores(B, edge_set, left_index_map, right_index_map, pool);

        for (const auto &l: *isolated_left_vertex_set) {
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }

        for (const auto &r: *isolated_right_vertex_set) {
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }
    }

    void branch_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                   const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                   const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                   const shared_ptr<thread_pool> &pool) {
        uint32_t max_k = 0;
        for(const auto &[k, k_list]:*core_order_index->get_middle_map()){
            if(!k_list->empty() && k > max_k){
                max_k = k;
            }
        }

        auto isolated_left_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_right_vertex_set = make_shared<unordered_set<uint32_t>>();

        B->remove_edge_collection(edge_set, isolated_left_vertex_set, isolated_right_vertex_set);

        auto global_left_mutex = make_shared<mutex>();
        auto global_right_mutex = make_shared<mutex>();


        for (uint32_t k = max_k; k >= 2; --k) {
            pool->submit_task([=] {
                /**
                 * @brief left path
                 */
                {
                    uint32_t max_i = k;
                    for (const auto &[i, i_list]: *core_order_index->get_left_map(k)) {
                        if (!i_list->empty() && i > max_i) {
                            max_i = i;
                        }
                    }

                    auto previous_removed_vector = make_shared<vector<uint32_t>>();
                    auto previous_removed_l_set = make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set = make_shared<unordered_set<uint32_t>>();

                    auto k_core_order_map = core_order_index->get_left_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_left_map(k);
                    auto k_core_degree_map = core_degree_index->get_left_map(k);

                    for (uint32_t i = max_i; i >= k; --i) {
                        auto removed_vector = make_shared<vector<uint32_t>>();

                        left_removal_partial_core(B, edge_set, k_core_order_map, k_core_rem_degree_map,
                                                  k_core_degree_map,
                                                  left_index_map, right_index_map, previous_removed_vector, i, k,
                                                  removed_vector);

                        for (const auto &v: *removed_vector) {
                            if (left_mutex_map->count(v)) {
                                left_mutex_map->at(v)->lock();
                                new_left_index_map->at(v)->remove(i, k);
                                left_mutex_map->at(v)->unlock();
                            } else {
                                right_mutex_map->at(v)->lock();
                                new_right_index_map->at(v)->remove(k, i);
                                right_mutex_map->at(v)->unlock();
                            }
                        }

                        swap(*previous_removed_vector, *removed_vector);
                    }

                    for (const auto &v: *previous_removed_vector) {
                        k_core_degree_map->erase(v);
                        k_core_rem_degree_map->erase(v);
                    }
                }
                /**
                 * @brief right path
                 */
                {
                    uint32_t max_j = k;
                    for (const auto &[j, j_list]: *core_order_index->get_right_map(k)) {
                        if (!j_list->empty() && j > max_j) {
                            max_j = j;
                        }
                    }

                    auto previous_removed_vector = make_shared<vector<uint32_t>>();

                    auto previous_removed_l_set = make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set = make_shared<unordered_set<uint32_t>>();

                    auto k_core_order_map = core_order_index->get_right_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_right_map(k);
                    auto k_core_degree_map = core_degree_index->get_right_map(k);

                    for (uint32_t j = max_j; j >= k; --j) {
                        auto removed_vector = make_shared<vector<uint32_t>>();

                        right_removal_partial_core(B, edge_set, k_core_order_map, k_core_rem_degree_map,
                                                   k_core_degree_map,
                                                   left_index_map, right_index_map, previous_removed_vector, k, j,
                                                   removed_vector);

                        for (const auto &v: *removed_vector) {
                            if (left_mutex_map->count(v)) {
                                left_mutex_map->at(v)->lock();
                                new_left_index_map->at(v)->remove(k, j);
                                left_mutex_map->at(v)->unlock();
                            } else {
                                right_mutex_map->at(v)->lock();
                                new_right_index_map->at(v)->remove(j, k);
                                right_mutex_map->at(v)->unlock();
                            }
                        }

                        swap(*previous_removed_vector, *removed_vector);
                    }

                    for (const auto &v: *previous_removed_vector) {
                        k_core_degree_map->erase(v);
                        k_core_rem_degree_map->erase(v);
                    }
                }
            });
        }
        {
            auto previous_middle_removed_vector = make_shared<vector<uint32_t>>();
            auto previous_middle_removed_l_set = make_shared<unordered_set<uint32_t>>();
            auto previous_middle_removed_r_set = make_shared<unordered_set<uint32_t>>();

            auto k_core_order_map = core_order_index->get_middle_map();
            auto k_core_rem_degree_map = core_rem_degree_index->get_middle_map();
            auto k_core_degree_map = core_degree_index->get_middle_map();

            for (uint32_t k = max_k; k >= 1; --k) {
                auto removed_vector = make_shared<vector<uint32_t>>();

                middle_removal_partial_core(B, edge_set, k_core_order_map, k_core_rem_degree_map, k_core_degree_map,
                                            left_index_map, right_index_map,
                                            previous_middle_removed_vector, k,
                                            removed_vector);

                for (const auto &v: *removed_vector) {
                    if (left_mutex_map->count(v)) {
                        left_mutex_map->at(v)->lock();
                        new_left_index_map->at(v)->remove(k, k);
                        left_mutex_map->at(v)->unlock();
                    } else {
                        right_mutex_map->at(v)->lock();
                        new_right_index_map->at(v)->remove(k, k);
                        right_mutex_map->at(v)->unlock();
                    }
                }

                swap(*previous_middle_removed_vector, *removed_vector);
            }

            for (const auto &v: *previous_middle_removed_vector) {
                k_core_degree_map->erase(v);
                k_core_rem_degree_map->erase(v);
            }
        }

        pool->barrier();


        for (const auto &[l, l_index]: *new_left_index_map) {
            pool->submit_task([=]{
                for (const auto &[i, j]: *l_index->get_index_map()) {
                    left_index_map->at(l)->remove(i, j - 1);
                    if (left_index_map->at(l)->get_j(i) == 0) {
                        left_index_map->at(l)->remove(i);
                    }
                }
                l_index->clear();
            });
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            pool->submit_task([=]{
                for (const auto [j, i]: *r_index->get_index_map()) {
                    right_index_map->at(r)->remove(j, i - 1);
                    if (right_index_map->at(r)->get_i(j) == 0) {
                        right_index_map->at(r)->remove(j);
                    }
                }
                r_index->clear();
            });
        }
        pool->barrier();

        update_trivial_bipartite_cores(B, edge_set, left_index_map, right_index_map, pool);

        for (const auto &l: *isolated_left_vertex_set) {
            left_index_map->erase(l);
            new_left_index_map->erase(l);
            left_mutex_map->erase(l);
        }

        for (const auto &r: *isolated_right_vertex_set) {
            right_index_map->erase(r);
            new_right_index_map->erase(r);
            right_mutex_map->erase(r);
        }
    }

    void branch_bipartite_core_maintenance::remove2(const shared_ptr<abstract_bipartite_graph> &B,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                    const shared_ptr<bipartite_core_order_index> &core_order_index,
                                                    const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
                                                    const shared_ptr<bipartite_core_degree_index> &core_degree_index,
                                                    const shared_ptr<thread_pool> &pool) {
        uint32_t max_k = 0;
        for (const auto &[k, k_list]: *core_order_index->get_middle_map()) {
            if (!k_list->empty() && k > max_k) {
                max_k = k;
            }
        }

        auto isolated_left_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_right_vertex_set = make_shared<unordered_set<uint32_t>>();

        B->remove_edge_collection(edge_set, isolated_left_vertex_set, isolated_right_vertex_set);

        auto global_left_mutex = make_shared<mutex>();
        auto global_right_mutex = make_shared<mutex>();


        for (uint32_t k = max_k; k >= 2; --k) {
            pool->submit_task([=] {
                /**
                 * @brief left path
                 */
                {
                    uint32_t max_i = k;
                    for (const auto &[i, i_list]: *core_order_index->get_left_map(k)) {
                        if (!i_list->empty() && i > max_i) {
                            max_i = i;
                        }
                    }

                    auto previous_removed_vector = make_shared<vector<uint32_t>>();
                    auto previous_removed_l_set = make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set = make_shared<unordered_set<uint32_t>>();

                    auto k_core_order_map = core_order_index->get_left_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_left_map(k);
                    auto k_core_degree_map = core_degree_index->get_left_map(k);

                    for (uint32_t i = max_i; i >= k; --i) {
                        auto removed_vector = make_shared<vector<uint32_t>>();

                        left_removal_partial_core(B, edge_set, k_core_order_map, k_core_rem_degree_map,
                                                  k_core_degree_map,
                                                  left_index_map, right_index_map, previous_removed_vector, i, k,
                                                  removed_vector);

                        for (const auto &v: *removed_vector) {
                            if (left_mutex_map->count(v)) {
                                left_mutex_map->at(v)->lock();
                                new_left_index_map->at(v)->remove(i, k);
                                left_mutex_map->at(v)->unlock();
                            } else {
                                right_mutex_map->at(v)->lock();
                                new_right_index_map->at(v)->remove(k, i);
                                right_mutex_map->at(v)->unlock();
                            }
                        }

                        swap(*previous_removed_vector, *removed_vector);
                    }

                    for (const auto &v: *previous_removed_vector) {
                        k_core_degree_map->erase(v);
                        k_core_rem_degree_map->erase(v);
                    }
                }
            });
        }

        for (uint32_t k = max_k; k >= 2; --k) {
            pool->submit_task([=] {
                /**
                 * @brief right path
                 */
                {
                    uint32_t max_j = k;
                    for (const auto &[j, j_list]: *core_order_index->get_right_map(k)) {
                        if (!j_list->empty() && j > max_j) {
                            max_j = j;
                        }
                    }

                    auto previous_removed_vector = make_shared<vector<uint32_t>>();

                    auto previous_removed_l_set = make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set = make_shared<unordered_set<uint32_t>>();

                    auto k_core_order_map = core_order_index->get_right_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_right_map(k);
                    auto k_core_degree_map = core_degree_index->get_right_map(k);

                    for (uint32_t j = max_j; j >= k; --j) {
                        auto removed_vector = make_shared<vector<uint32_t>>();

                        right_removal_partial_core(B, edge_set, k_core_order_map, k_core_rem_degree_map,
                                                   k_core_degree_map,
                                                   left_index_map, right_index_map, previous_removed_vector, k, j,
                                                   removed_vector);

                        for (const auto &v: *removed_vector) {
                            if (left_mutex_map->count(v)) {
                                left_mutex_map->at(v)->lock();
                                new_left_index_map->at(v)->remove(k, j);
                                left_mutex_map->at(v)->unlock();
                            } else {
                                right_mutex_map->at(v)->lock();
                                new_right_index_map->at(v)->remove(j, k);
                                right_mutex_map->at(v)->unlock();
                            }
                        }

                        swap(*previous_removed_vector, *removed_vector);
                    }

                    for (const auto &v: *previous_removed_vector) {
                        k_core_degree_map->erase(v);
                        k_core_rem_degree_map->erase(v);
                    }
                }
            });
        }
        {
            auto previous_middle_removed_vector = make_shared<vector<uint32_t>>();
            auto previous_middle_removed_l_set = make_shared<unordered_set<uint32_t>>();
            auto previous_middle_removed_r_set = make_shared<unordered_set<uint32_t>>();

            auto k_core_order_map = core_order_index->get_middle_map();
            auto k_core_rem_degree_map = core_rem_degree_index->get_middle_map();
            auto k_core_degree_map = core_degree_index->get_middle_map();

            for (uint32_t k = max_k; k >= 1; --k) {
                auto removed_vector = make_shared<vector<uint32_t>>();

                middle_removal_partial_core(B, edge_set, k_core_order_map, k_core_rem_degree_map, k_core_degree_map,
                                            left_index_map, right_index_map,
                                            previous_middle_removed_vector, k,
                                            removed_vector);

                for (const auto &v: *removed_vector) {
                    if (left_mutex_map->count(v)) {
                        left_mutex_map->at(v)->lock();
                        new_left_index_map->at(v)->remove(k, k);
                        left_mutex_map->at(v)->unlock();
                    } else {
                        right_mutex_map->at(v)->lock();
                        new_right_index_map->at(v)->remove(k, k);
                        right_mutex_map->at(v)->unlock();
                    }
                }

                swap(*previous_middle_removed_vector, *removed_vector);
            }

            for (const auto &v: *previous_middle_removed_vector) {
                k_core_degree_map->erase(v);
                k_core_rem_degree_map->erase(v);
            }
        }

        pool->barrier();


        for (const auto &[l, l_index]: *new_left_index_map) {
            pool->submit_task([=] {
                for (const auto &[i, j]: *l_index->get_index_map()) {
                    left_index_map->at(l)->remove(i, j - 1);
                    if (left_index_map->at(l)->get_j(i) == 0) {
                        left_index_map->at(l)->remove(i);
                    }
                }
                l_index->clear();
            });
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            pool->submit_task([=] {
                for (const auto [j, i]: *r_index->get_index_map()) {
                    right_index_map->at(r)->remove(j, i - 1);
                    if (right_index_map->at(r)->get_i(j) == 0) {
                        right_index_map->at(r)->remove(j);
                    }
                }
                r_index->clear();
            });
        }
        pool->barrier();

        update_trivial_bipartite_cores(B, edge_set, left_index_map, right_index_map, pool);

        for (const auto &l: *isolated_left_vertex_set) {
            left_index_map->erase(l);
            new_left_index_map->erase(l);
            left_mutex_map->erase(l);
        }

        for (const auto &r: *isolated_right_vertex_set) {
            right_index_map->erase(r);
            new_right_index_map->erase(r);
            right_mutex_map->erase(r);
        }
    }

    void branch_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                   const shared_ptr<thread_pool> &pool) {
        uint32_t max_k = 0;
        {
            while(true){
                auto flag = false;
                for (const auto &e: *edge_set) {
                    auto l = e->get_left_vertex_id();
                    auto r = e->get_right_vertex_id();

                    if(left_index_map->at(l)->count(max_k + 1, max_k + 1) && right_index_map->at(r)->count(max_k + 1, max_k + 1)){
                        flag = true;
                        break;
                    }
                }
                if(flag){
                    ++max_k;
                }else{
                    break;
                }
            }

            for(const auto &e:*edge_set){
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                auto l_degree = B->get_left_vertex(l)->get_degree();
                right_core_degree_map->at(r)->at(l_degree)->erase(l);
                if(right_core_degree_map->at(r)->empty()){
                    right_core_degree_map->at(r)->erase(l_degree);
                }

                auto r_degree = B->get_right_vertex(r)->get_degree();
                left_core_degree_map->at(l)->at(r_degree)->erase(r);
                if (left_core_degree_map->at(l)->at(r_degree)->empty()) {
                    left_core_degree_map->at(l)->erase(r_degree);
                }
            }
        }

        auto isolated_l_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_r_set = make_shared<unordered_set<uint32_t>>();
        B->remove_edge_collection(edge_set, isolated_l_set, isolated_r_set);

        auto global_left_mutex = make_shared<mutex>();
        auto global_right_mutex = make_shared<mutex>();

        for(uint32_t k = max_k; k >= 2; --k){
            {
                pool->submit_task([=]{
                    uint32_t max_i = k;
                    while (true) {
                        auto flag = false;
                        for (const auto &e: *edge_set) {
                            auto l = e->get_left_vertex_id();
                            auto r = e->get_right_vertex_id();

                            if (left_index_map->at(l)->count(max_i + 1, k) &&
                                right_index_map->at(r)->count(k, max_i + 1)) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) {
                            ++max_i;
                        } else {
                            break;
                        }
                    }

                    auto previous_removed_l_set =  make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set =  make_shared<unordered_set<uint32_t>>();

                    for (uint32_t  i = max_i; i > k; --i){
                        auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                        auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                        left_removal_partial_core(B, left_core_degree_map, right_core_degree_map, edge_set, left_index_map, right_index_map,
                                                  previous_removed_l_set, previous_removed_r_set, i, k, removed_l_set, removed_r_set);

                        global_left_mutex->lock();
                        for(const auto&l:*removed_l_set){
                            new_left_index_map->at(l)->remove(i, k);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for(const auto&r:*removed_r_set){
                            new_right_index_map->at(r)->remove(k, i);
                        }
                        global_right_mutex->unlock();

                        swap(*previous_removed_l_set, *removed_l_set);
                        swap(*previous_removed_r_set, *removed_r_set);
                    }
                });
            }

            {
                pool->submit_task([=]{
                    uint32_t max_j = k;
                    while(true){
                        auto flag = false;
                        for (const auto &e: *edge_set) {
                            auto l = e->get_left_vertex_id();
                            auto r = e->get_right_vertex_id();

                            if (left_index_map->at(l)->count(k, max_j + 1) &&
                                right_index_map->at(r)->count(max_j + 1, k)) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) {
                            ++max_j;
                        } else {
                            break;
                        }
                    }

                    auto previous_removed_l_set =  make_shared<unordered_set<uint32_t>>();
                    auto previous_removed_r_set =  make_shared<unordered_set<uint32_t>>();

                    for (uint32_t  j = max_j; j > k; --j){
                        auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                        auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                        right_removal_partial_core(B, left_core_degree_map, right_core_degree_map, edge_set, left_index_map, right_index_map,
                                                   previous_removed_l_set, previous_removed_r_set,k, j, removed_l_set, removed_r_set);
                        global_left_mutex->lock();
                        for(const auto &l:*removed_l_set){
                            new_left_index_map->at(l)->remove(k, j);
                        }
                        global_left_mutex->unlock();

                        global_right_mutex->lock();
                        for(const auto &r:*removed_r_set){
                            new_right_index_map->at(r)->remove(j, k);
                        }
                        global_right_mutex->unlock();

                        swap(*previous_removed_l_set, *removed_l_set);
                        swap(*previous_removed_r_set, *removed_r_set);
                    }
                });
            }
        }

        auto previous_middle_removed_l_set =  make_shared<unordered_set<uint32_t>>();
        auto previous_middle_removed_r_set =  make_shared<unordered_set<uint32_t>>();

        for(uint32_t k = max_k; k >= 2; --k){
            auto removed_l_set = make_shared<unordered_set<uint32_t>>();
            auto removed_r_set = make_shared<unordered_set<uint32_t>>();

            middle_removal_partial_core(B, left_core_degree_map, right_core_degree_map, edge_set, left_index_map, right_index_map,
                                        previous_middle_removed_l_set, previous_middle_removed_r_set, k, removed_l_set, removed_r_set);

            global_left_mutex->lock();
            for(const auto &l:*removed_l_set){
                new_left_index_map->at(l)->remove(k, k);
            }
            global_left_mutex->unlock();

            global_right_mutex->lock();
            for(const auto &r:*removed_r_set){
                new_right_index_map->at(r)->remove(k, k);
            }
            global_right_mutex->unlock();

            swap(*previous_middle_removed_l_set, *removed_l_set);
            swap(*previous_middle_removed_r_set, *removed_r_set);
        }

        pool->barrier();

        update_index_for_removal(B, edge_set, left_core_degree_map, right_core_degree_map,  isolated_l_set, isolated_r_set);

        /**
         * @brief update core index
         */
        for(const auto&[l, l_index]:*new_left_index_map) {
            for (const auto &[i, j]: *l_index->get_index_map()) {
                left_index_map->at(l)->remove(i, j - 1);
                if(left_index_map->at(l)->get_j(i) == 0){
                    left_index_map->at(l)->remove(i);
                }
            }
            l_index->clear();
        }

        for(const auto &[r, r_index]:*new_right_index_map){
            for(const auto [j, i]:*r_index->get_index_map()){
                right_index_map->at(r)->remove(j, i - 1);
                if( right_index_map->at(r)->get_i(j) == 0){
                    right_index_map->at(r)->remove(j);
                }
            }
            r_index->clear();
        }

        update_trivial_bipartite_cores(B, edge_set, left_index_map, right_index_map, pool);


        for(const auto&l:*isolated_l_set){
            left_core_degree_map->erase(l);
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }

        for(const auto &r:*isolated_r_set){
            right_core_degree_map->erase(r);
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }
    }


    void branch_bipartite_core_maintenance::candidate_graph(const shared_ptr<abstract_bipartite_graph>& B,
                                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                            const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                            uint32_t k,
                                                            const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                                                            const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map){
        auto visited_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto visited_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set  = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief get (i,j)-insert graph
         */
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            /**
             * @brief l is not in the previous (i,j)-core
             */
            if (!left_index_map->at(l)->count(k, k)) {
                auto core_degree = B->get_left_vertex(l)->get_degree();
                if (core_degree >= k) {
                    if(visited_l_map->count(l)){
                        ++visited_l_map->at(l);
                    }else{
                        visited_l_map->insert({l, 1});
                    }
                } else {
                    evicted_l_set->insert(l);
                }
            }

            /**
             * @brief r is not in the previous (i,j)-core
             */
            if (!right_index_map->at(r)->count(k, k)) {
                auto core_degree = B->get_right_vertex(r)->get_degree();

                if (core_degree >= k) {
                    if(visited_r_map->count(r)){
                        ++visited_r_map->at(r);
                    }else{
                        visited_r_map->insert({r, 1});
                    }
                } else {
                    evicted_r_set->insert(r);
                }
            }
        }

        while(!visited_l_map->empty() || !visited_r_map->empty()) {
            while (!visited_l_map->empty()) {
                auto [l, l_count] = *visited_l_map->begin();
                visited_l_map->erase(l);

                if (l_count == 0 || evicted_l_set->count(l)) {
                    continue;
                }

                uint32_t degree = 0;
                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (inserted_r_map->count(r) || right_index_map->at(r)->count(k, k)) {
                        ++degree;
                    } else if (!evicted_r_set->count(r)) {
                        auto core_degree = B->get_right_vertex(r)->get_degree();

                        if (core_degree >= k) {
                            ++degree;
                            if(visited_r_map->count(r)){
                                ++visited_r_map->at(r);
                            }else{
                                visited_r_map->insert({r, 1});
                            }
                        } else {
                            evicted_r_set->insert(r);
                        }
                    }
                }

                if (degree >= k) {
                    inserted_l_map->insert({l, degree});
                } else {
                    remove_left_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                       visited_l_map, visited_r_map, l, k, k);
                }
            }

            while (!visited_r_map->empty()) {
                auto [r, r_count] = *visited_r_map->begin();
                visited_r_map->erase(r);

                if (r_count == 0 || evicted_r_set->count(r)) {
                    continue;
                }

                uint32_t degree = 0;
                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (inserted_l_map->count(l) || left_index_map->at(l)->count(k, k)) {
                        ++degree;
                    } else if (!evicted_l_set->count(l)) {
                        auto core_degree = B->get_left_vertex(l)->get_degree();

                        if (core_degree >= k) {
                            ++degree;
                            if(visited_l_map->count(l)){

                                ++visited_l_map->at(l);
                            }else{
                                visited_l_map->insert({l, 1});
                            }
                        } else {
                            evicted_l_set->insert(l);
                        }
                    }
                }
                if (degree >= k) {
                    inserted_r_map->insert({r, degree});
                } else {
                    remove_right_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                        visited_l_map, visited_r_map, r, k, k);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_order_map,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_degree_map,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                            uint32_t k,
                                                            const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                            const shared_ptr<unordered_set<uint32_t>> &inserted_r_set) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        /**
         * @brief get (i,j)-insert graph
         */
        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(!left_index_map->at(l)->count(k - 1, k - 1) || !right_index_map->at(r)->count(k - 1, k - 1)){
                continue;
            }

            if(left_index_map->at(l)->count(k - 1, k - 1) && !left_index_map->at(l)->count(k, k)){
                inserted_l_set->insert(l);
            }

            if(right_index_map->at(r)->count(k - 1, k - 1) && !right_index_map->at(r)->count(k, k) ){
                inserted_r_set->insert(r);
            }
        }

        for(const auto &l:*inserted_l_set){
            vertex_degree_map->insert({l, B->get_left_vertex(l)->get_degree()});
        }

        for(const auto &r:*inserted_r_set){
            vertex_degree_map->insert({r, B->get_right_vertex(r)->get_degree()});
        }

        auto insertion_vertex_vector = make_shared<vector<uint32_t>>();
        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                insertion_vertex_vector->push_back(l);

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if(vertex_degree_map->count(r)){
                        --vertex_degree_map->at(r);
                        if(vertex_degree_map->at(r) == k - 1){
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                insertion_vertex_vector->push_back(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (vertex_degree_map->count(l)) {
                        --vertex_degree_map->at(l);
                        if(vertex_degree_map->at(l) == k - 1){
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }

        auto order_list = k_order_map->at(k - 1);
        auto affected_vertex_map = make_shared<map<uint32_t, uint32_t>>();
        auto ext = make_shared<unordered_map<uint32_t, uint32_t>>();

        for(const auto &l:*inserted_l_set){
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if(order_list->count_value(r)){
                    if(!ext->count(r)){
                        ext->insert({r, 0});
                        affected_vertex_map->insert({order_list->find_key(r).value(), r});
                    }
                    ++ext->at(r);
                }
            }
        }

        for(const auto &r:*inserted_r_set){
            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if(order_list->count_value(l)){
                    if(!ext->count(l)){
                        ext->insert({l, 0});
                        affected_vertex_map->insert({order_list->find_key(l).value(), l});
                    }
                    ++ext->at(l);
                }
            }
        }

        for(const auto &v:*insertion_vertex_vector){
            order_list->left_insert(v);
        }

        while(!affected_vertex_map->empty()){
            auto [key, v] = *affected_vertex_map->begin();
            affected_vertex_map->erase(key);

            if(ext->at(v) + k_rem_degree_map->at(v) >= k){
                vertex_degree_map->insert({v, ext->at(v) + k_rem_degree_map->at(v)});
                ext->at(v) = 0;

                order_list->remove(v);

                if(B->get_left_vertex(v)){
                    auto &l = v;
                    inserted_l_set->insert(l);
                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if(order_list->count_value(r) && order_list->find_key(r).value() > key){
                            if(!ext->count(r)){
                                ext->insert({r, 0});
                                affected_vertex_map->insert({order_list->find_key(r).value(), r});
                            }
                            ++ext->at(r);
                        }
                    }
                }else{
                    auto &r = v;
                    inserted_r_set->insert(r);
                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if(order_list->count_value(l) && order_list->find_key(l).value() > key){
                            if(!ext->count(l)){
                                ext->insert({l, 0});
                                affected_vertex_map->insert({order_list->find_key(l).value(), l});
                            }
                            ++ext->at(l);
                        }
                    }
                }
            }else{
                k_rem_degree_map->at(v) = ext->at(v) + k_rem_degree_map->at(v);
                ext->at(v) = 0;

                if(B->get_left_vertex(v)){
                    auto &l = v;
                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if(vertex_degree_map->count(r)){
                            --vertex_degree_map->at(r);
                            if(vertex_degree_map->at(r) < k){
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }else{
                    auto &r = v;
                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if(vertex_degree_map->count(l)){
                            --vertex_degree_map->at(l);
                            if(vertex_degree_map->at(l) < k){
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }

                auto v_next_node = order_list->find(v)->get_next();
                while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                    while (!evicted_l_set->empty()) {
                        auto l = *evicted_l_set->begin();
                        evicted_l_set->erase(l);

                        auto l_node = make_shared<extend_node<int, uint32_t>>(0, l);
                        order_list->insert_before(l_node, v_next_node);
                        k_rem_degree_map->at(l) = vertex_degree_map->at(l);

                        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                            if(vertex_degree_map->count(r)){
                                --vertex_degree_map->at(r);
                                if(vertex_degree_map->at(r) < k){
                                    evicted_r_set->insert(r);
                                }
                            }else if(ext->count(r)){
                                --ext->at(r);
                                if(ext->at(r) == 0){
                                    ext->erase(r);
                                    affected_vertex_map->erase(order_list->find_key(r).value());
                                }
                            }
                        }
                    }

                    while (!evicted_r_set->empty()) {
                        auto r = *evicted_r_set->begin();
                        evicted_r_set->erase(r);

                        auto r_node = make_shared<extend_node<int, uint32_t>>(0, r);
                        order_list->insert_before(r_node, v_next_node);
                        k_rem_degree_map->at(r) = vertex_degree_map->at(r);

                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (vertex_degree_map->count(l)) {
                                --vertex_degree_map->at(l);
                                if(vertex_degree_map->at(l) < k){
                                    evicted_l_set->insert(l);
                                }
                            }else if(ext->count(l)){
                                --ext->at(l);
                                if(ext->at(l) == 0){
                                    ext->erase(l);
                                    affected_vertex_map->erase(order_list->find_key(l).value());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    uint32_t
    branch_bipartite_core_maintenance::compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &B,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_vertex_index_map,
                                                                       uint32_t l,
                                                                       uint32_t i,
                                                                       uint32_t j) {
        uint32_t core_degree = 0;
        if(!B->get_left_vertex(l)){
            return core_degree;
        }
        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
            if(right_vertex_index_map->at(r)->count(j, i)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph>&B,
                                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_vertex_index_map,
                                                                                 uint32_t r,
                                                                                 uint32_t i,
                                                                                 uint32_t j){
        uint32_t core_degree = 0;
        if(!B->get_right_vertex(r)){
            return 0;
        }

        for(const auto&[l,e]:*B->get_right_vertex(r)->get_edge_map()){
            if(left_vertex_index_map->at(l)->count(i, j)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_left_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                                                uint32_t l,
                                                                                uint32_t i,
                                                                                uint32_t j)
    {
        uint32_t core_degree = 0;
        auto degree_map = left_core_degree_map->at(l);
        for(auto iter = degree_map->lower_bound(j); iter!=degree_map->end();++iter){
            for(const auto&r:*iter->second){
                if(right_index_map->at(r)->count(j, i)){
                    ++core_degree;
                }
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_right_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                                 uint32_t r,
                                                                                 uint32_t i,
                                                                                 uint32_t j){
        uint32_t core_degree = 0;
        auto degree_map = right_core_degree_map->at(r);
        for(auto iter = degree_map->lower_bound(i); iter !=degree_map->end(); ++iter){
            for(const auto&l:*iter->second){
                if(left_index_map->at(l)->count(i, j)){
                    ++core_degree;
                }
            }
        }

        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph>&B,
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_vertex_index_map,
                                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                                                uint32_t l,
                                                                                uint32_t i,
                                                                                uint32_t j)
    {
        uint32_t core_degree = 0;
        for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
            if(inserted_r_map->count(r) || right_vertex_index_map->at(r)->count(j, i)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph>&B,
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_vertex_index_map,
                                                                                const shared_ptr<unordered_set<uint32_t>>& removed_r_set,
                                                                                uint32_t l,
                                                                                uint32_t i,
                                                                                uint32_t j)
    {
        uint32_t core_degree = 0;
        for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
            if(!removed_r_set->count(r) &&  right_vertex_index_map->at(r)->count(j, i)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_left_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_vertex_index_map,
                                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                                                uint32_t l,
                                                                                uint32_t i,
                                                                                uint32_t j)
    {
        uint32_t core_degree = 0;
        auto degree_map = left_core_degree_map->at(l);
        for(auto iter = degree_map->lower_bound(j); iter != degree_map->end(); ++iter){
            for(const auto &r:*iter->second){
                if(inserted_r_map->count(r) || right_vertex_index_map->at(r)->count(j, i)){
                    ++core_degree;
                }
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph>&B,
                                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_vertex_index_map,
                                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                                                 uint32_t r,
                                                                                 uint32_t i,
                                                                                 uint32_t j){
        uint32_t core_degree = 0;
        for(const auto&[l,e]:*B->get_right_vertex(r)->get_edge_map()){
            if(inserted_l_map->count(l) || left_vertex_index_map->at(l)->count(i, j)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph>&B,
                                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_vertex_index_map,
                                                                                 const shared_ptr<unordered_set<uint32_t>>& removed_l_set,
                                                                                 uint32_t r,
                                                                                 uint32_t i,
                                                                                 uint32_t j){
        uint32_t core_degree = 0;
        for(const auto&[l,e]:*B->get_right_vertex(r)->get_edge_map()){
            if(!removed_l_set->count(l) && left_vertex_index_map->at(l)->count(i, j)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_right_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>>& left_vertex_index_map,
                                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                                                 uint32_t r,
                                                                                 uint32_t i,
                                                                                 uint32_t j){
        uint32_t core_degree = 0;
        auto degree_map = right_core_degree_map->at(r);
        for(auto iter = degree_map->lower_bound(i); iter!=degree_map->end(); ++iter){
            for(const auto &l:*iter->second){
                if(inserted_l_map->count(l) || left_vertex_index_map->at(l)->count(i, j)){
                    ++core_degree;
                }
            }
        }
        return core_degree;
    }

    void branch_bipartite_core_maintenance::left_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                                                                 uint32_t i,
                                                                 uint32_t k,
                                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map){
        auto visited_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto visited_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set  = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        /**
         * @brief get (i,k)-insert graph
         */
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            /**
             * @brief skip edges are not in the current (i-1,k)-core
             */
            if((previous_inserted_l_map->count(l) || left_index_map->at(l)->count(i - 1 , k))
               && (previous_inserted_r_map->count(r) || right_index_map->at(r)->count(k, i - 1))){
                /**
                  * @brief l is not in the previous (i,k)-core
                  */
                if(!left_index_map->at(l)->count(i, k))
                {
                    if(!l_map->count(l)){
                        auto core_degree = previous_inserted_l_map->count(l) ? previous_inserted_l_map->at(l):compute_left_vertex_core_degree(B, right_index_map, previous_inserted_r_map, l, i - 1, k);
                        l_map->insert({l, core_degree});
                    }
                    if(l_map->at(l) >= i){
                        if(visited_l_map->count(l)){
                            ++visited_l_map->at(l);
                        }else{
                            visited_l_map->insert({l,1});
                        }
                    }else{
                        evicted_l_set->insert(l);
                    }
                }

                /**
                 * @brief r is not in the previous (i,k)-core
                 */
                if(!right_index_map->at(r)->count(k, i)){
                    if(visited_r_map->count(r)) {
                        ++visited_r_map->at(r);
                    }else{
                        visited_r_map->insert({r, 1});
                    }
                }
            }
        }

        while(!visited_l_map->empty() || !visited_r_map->empty()){
            while(!visited_l_map->empty()){
                auto [l, l_count] = *visited_l_map->begin();
                visited_l_map->erase(l);

                if(l_count == 0 || evicted_l_set->count(l)){
                    continue;
                }

                uint32_t  degree = 0;
                for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if (!(right_index_map->at(r)->count(k, i-1) || previous_inserted_r_map->count(r))) {
                        continue;
                    }
                    if (right_index_map->at(r)->count(k, i) || inserted_r_map->count(r)) {
                        ++degree;
                    } else if (!evicted_r_set->count(r)) {
                        ++degree;
                        if(visited_r_map->count(r)){
                            ++visited_r_map->at(r);
                        }else{
                            visited_r_map->insert({r, 1});
                        }
                    }
                }
                if(degree >= i){
                    inserted_l_map->insert({l, degree});
                } else
                {
                    remove_left_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                       visited_l_map, visited_r_map,l, i, k);
                }
            }

            while(!visited_r_map->empty()){
                auto [r, r_count] = *visited_r_map->begin();
                visited_r_map->erase(r);

                if(r_count == 0 || evicted_r_set->count(r)){
                    continue;
                }

                uint32_t degree = 0;
                for(const auto& [l, e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(!(left_index_map->at(l)->count(i - 1, k)|| previous_inserted_l_map->count(l))){
                        continue;
                    }
                    if(left_index_map->at(l)->count(i, k) || inserted_l_map->count(l)){
                        ++degree;
                    }else  if(!evicted_l_set->count(l)){
                        if(!l_map->count(l)){
                            auto core_degree = previous_inserted_l_map->count(l) ?
                                               previous_inserted_l_map->at(l): compute_left_vertex_core_degree(B, right_index_map, previous_inserted_r_map, l, i - 1, k);
                            l_map->insert({l, core_degree});
                        }

                        if(l_map->at(l) >= i){
                            ++degree;
                            if(visited_l_map->count(l)){
                                ++visited_l_map->at(l);
                            }else{
                                visited_l_map->insert({l, 1});
                            }
                        } else{
                            evicted_l_set->insert(l);
                        }
                    }
                }
                if(degree >= k){
                    inserted_r_map->insert({r, degree});
                } else {
                    remove_right_vertex(B, inserted_l_map, inserted_r_map,
                                        evicted_l_set, evicted_r_set, visited_l_map, visited_r_map, r, i, k);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::left_candidate_graph(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                 const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                                 const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_order_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                 uint32_t i,
                                                                 uint32_t k) {
        if(!k_order_map->count(i - 1)){
            k_order_map->insert({i - 1, make_shared<extend_list<int, uint32_t>>()});
        }
        auto order_list = k_order_map->at(i - 1);

//        for(auto p = order_list->get_head(); p ; p = p->get_next()){
//            printf("%u ", p->get_value());
//        }
//        printf("\n");

        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if((!left_index_map->at(l)->count(i - 1, k) && !inserted_l_set->count(l)) ||
               (!right_index_map->at(r)->count(k, i - 1) && !inserted_r_set->count(r)))
            {
                continue;
            }

            if (left_index_map->at(l)->count(i, k) && right_index_map->at(r)->count(k, i)
                && !(left_index_map->at(l)->count(i + 1, k) && right_index_map->at(r)->count(k, i + 1))) {

                if (left_index_map->at(l)->count(i, k) && !left_index_map->at(l)->count(i + 1, k)) {
                    ++k_core_degree_map->at(l);
                }

                if (right_index_map->at(r)->count(k, i) && !right_index_map->at(r)->count(k, i + 1)) {
                    ++k_core_degree_map->at(r);
                }
                continue;
            }

            if (order_list->count_value(l)) {
                auto key = order_list->find_key(l).value();
                for (const auto &[r1, e1]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (order_list->count_value(r1) && order_list->find_key(r1).value() < key) {
                        --k_rem_map->at(r1);
                    }
                }

                order_list->remove(l);
                inserted_l_set->insert(l);
            }

            if(order_list->count_value(r)){
                auto key = order_list->find_key(r).value();
                for (const auto &[l1, e1]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (order_list->count_value(l1) && order_list->find_key(l1).value() < key) {
                        --k_rem_map->at(l1);
                    }
                }

                order_list->remove(r);
                inserted_r_set->insert(r);
            }
        }

        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &l: *inserted_l_set) {
            vertex_degree_map->insert({l, k_core_degree_map->at(l)});
            if (k_core_degree_map->at(l) == i - 1) {
                evicted_l_set->insert(l);
            }
        }

        for (const auto &r: *inserted_r_set) {
            vertex_degree_map->insert({r, k_core_degree_map->at(r)});
        }

        auto insertion_vertex_vector = make_shared<vector<uint32_t>>();

        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                insertion_vertex_vector->push_back(l);
                inserted_l_set->erase(l);
                vertex_degree_map->erase(l);

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (vertex_degree_map->count(r)) {
                        --vertex_degree_map->at(r);
                        if (vertex_degree_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                insertion_vertex_vector->push_back(r);
                inserted_r_set->erase(r);
                vertex_degree_map->erase(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (vertex_degree_map->count(l)) {
                        --vertex_degree_map->at(l);
                        if (vertex_degree_map->at(l) == i - 1) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }

        auto ext = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto affected_vertex_map = make_shared<map<int, uint32_t>>();

        for (const auto &l: *inserted_l_set) {
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (order_list->count_value(r)) {
                    if (!ext->count(r)) {
                        ext->insert({r, 0});
                        affected_vertex_map->insert({order_list->find_key(r).value(), r});
                    }
                    ++ext->at(r);
                }
            }
        }

        for (const auto &r: *inserted_r_set) {
            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if (order_list->count_value(l)) {
                    if (!ext->count(l)) {
                        ext->insert({l, 0});
                        affected_vertex_map->insert({order_list->find_key(l).value(), l});
                    }
                    ++ext->at(l);
                }
            }
        }

//        for(auto p = order_list->get_head(); p; p = p->get_next()){
//            printf("%u ", p->get_value());
//        }
//        printf("\n");

        for(auto iter = insertion_vertex_vector->rbegin(); iter!=insertion_vertex_vector->rend(); ++iter){
            order_list->left_insert(*iter);
        }

        while (!affected_vertex_map->empty()) {
            auto [key, v] = *affected_vertex_map->begin();
            affected_vertex_map->erase(key);

            if (B->get_left_vertex(v)) {
                auto &l = v;
                if (ext->at(l) + k_rem_map->at(l) >= i) {
                    order_list->remove(l);
                    inserted_l_set->insert(l);
                    vertex_degree_map->insert({l, ext->at(l) + k_rem_map->at(l)});
                    ext->at(l) = 0;

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (order_list->count_value(r) && order_list->find_key(r).value() > key) {
                            if (!ext->count(r)) {
                                ext->insert({r, 0});
                                affected_vertex_map->insert({order_list->find_key(r).value(), r});
                            }
                            ++ext->at(r);
                        }
                    }
                } else {
                    k_rem_map->at(l) = ext->at(l) + k_rem_map->at(l);
                    ext->at(l) = 0;

                    for(const auto &[r, e]:*B->get_left_vertex(l)->get_edge_map()){
                        if(vertex_degree_map->count(r)){
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) < k) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }

                    remove_vertex(B, l, evicted_l_set, evicted_r_set, inserted_l_set, inserted_r_set, vertex_degree_map, order_list,k_rem_map, ext, affected_vertex_map, i, k);
                }
            } else {
                auto &r = v;
                if(ext->at(r) + k_rem_map->at(r) >= k){
                    order_list->remove(r);
                    inserted_r_set->insert(r);
                    vertex_degree_map->insert({r, ext->at(r) + k_rem_map->at(r)});
                    ext->at(r) = 0;

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (order_list->count_value(l) && order_list->find_key(l).value() > key) {
                            if (!ext->count(l)) {
                                ext->insert({l, 0});
                                affected_vertex_map->insert({order_list->find_key(l).value(), l});
                            }
                            ++ext->at(l);
                        }
                    }
                }else{
                    k_rem_map->at(r) = ext->at(r) + k_rem_map->at(r);
                    ext->at(r) = 0;

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (vertex_degree_map->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) == i - 1) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }

                    remove_vertex(B, v, evicted_l_set, evicted_r_set, inserted_l_set, inserted_r_set, vertex_degree_map,
                                  order_list, k_rem_map, ext, affected_vertex_map, i, k);
                }
            }
        }

        order_list->reset_order();

        /**
         * @brief update the core degree of updated vertices
         */
        for (const auto &l: *inserted_l_set) {
            k_core_degree_map->at(l) = vertex_degree_map->at(l);
        }

        for (const auto &r: *inserted_r_set) {
            k_core_degree_map->at(r) = vertex_degree_map->at(r);
        }
    }

    void branch_bipartite_core_maintenance::left_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
            uint32_t i,
            uint32_t k,
            const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
            const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map){
        auto visited_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto visited_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set  = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        /**
         * @brief get (i,k)-insert graph
         */
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            /**
             * @brief skip edges are not in the current (i-1,k)-core
             */
            if((previous_inserted_l_map->count(l) || left_index_map->at(l)->count(i - 1 , k))
               && (previous_inserted_r_map->count(r) || right_index_map->at(r)->count(k, i - 1))){
                /**
                 * @brief l is not in the previous (i,k)-core
                 */
                if(!left_index_map->at(l)->count(i, k))
                {
                    if(!l_map->count(l)){
                        auto core_degree = previous_inserted_l_map->count(l) ?
                                           previous_inserted_l_map->at(l):compute_left_vertex_core_degree(left_core_degree_map, right_index_map, previous_inserted_r_map, l, i - 1, k);
                        l_map->insert({l, core_degree});
                    }
                    if(l_map->at(l) >= i){
                        if(visited_l_map->count(l)){
                            ++visited_l_map->at(l);
                        } else{
                            visited_l_map->insert({l, 1});
                        }
                    }else{
                        evicted_l_set->count(l);
                    }
                }

                /**
                 * @brief r is not in the previous (i,k)-core
                 */
                if(!right_index_map->at(r)->count(k, i)){
                    if(visited_r_map->count(r)){
                        ++visited_r_map->at(r);
                    } else{
                        visited_r_map->insert({r, 1});
                    }
                }
            }
        }

        while(!visited_l_map->empty() || !visited_r_map->empty()){
            while(!visited_l_map->empty()){
                auto [l, l_count] = *visited_l_map->begin();
                visited_l_map->erase(l);

                if(l_count == 0 || evicted_l_set->count(l)){
                    continue;
                }
                uint32_t  degree = 0;
                auto degree_map = left_core_degree_map->at(l);
                for(auto iter = degree_map->lower_bound(k);iter!=degree_map->end();++iter){
                    for(const auto &r:*iter->second) {
                        if (!(right_index_map->at(r)->count(k, i-1) || previous_inserted_r_map->count(r))) {
                            continue;
                        }
                        if (right_index_map->at(r)->count(k, i) || inserted_r_map->count(r)) {
                            ++degree;
                        } else if (!evicted_r_set->count(r)) {
                            ++degree;
                            if(visited_r_map->count(r)){
                                ++visited_r_map->at(r);
                            }else{
                                visited_r_map->insert({r, 1});
                            }
                        }
                    }
                }

                if(degree >= i){
                    inserted_l_map->insert({l, degree});
                } else
                {
                    remove_left_vertex(left_core_degree_map, right_core_degree_map,
                                       inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                       visited_l_map, visited_r_map,
                                       l, i, k);
                }
            }



            while(!visited_r_map->empty()){
                auto [r, r_count] = *visited_r_map->begin();
                visited_r_map->erase(r);

                if(r_count == 0 || evicted_r_set->count(r)){
                    continue;
                }

                uint32_t degree = 0;
                auto degree_map = right_core_degree_map->at(r);
                for(auto iter = degree_map->lower_bound(i); iter!=degree_map->end();++iter){
                    for(const auto &l:*iter->second){
                        if(!(left_index_map->at(l)->count(i - 1, k)|| previous_inserted_l_map->count(l))){
                            continue;
                        }
                        if(left_index_map->at(l)->count(i, k) || inserted_l_map->count(l)){
                            ++degree;
                        }else  if(!evicted_l_set->count(l)){
                            if(!l_map->count(l)){
                                auto core_degree = previous_inserted_l_map->count(l) ?
                                                   previous_inserted_l_map->at(l): compute_left_vertex_core_degree(left_core_degree_map, right_index_map, previous_inserted_r_map, l, i - 1, k);
                                l_map->insert({l, core_degree});
                            }

                            if(l_map->at(l) >= i){
                                ++degree;
                                if(visited_l_map->count(l)){
                                    ++visited_l_map->at(l);
                                }else{
                                    visited_l_map->insert({l, 1});
                                }
                            }else{
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
                if(degree >= k){
                    inserted_r_map->insert({r, degree});
                } else
                {
                    remove_right_vertex(left_core_degree_map, right_core_degree_map, inserted_l_map, inserted_r_map,
                                        evicted_l_set, evicted_r_set,
                                        visited_l_map, visited_r_map,
                                        r, i, k);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::middle_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                                                                   uint32_t k,
                                                                   const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                                                                   const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map){
        auto l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto visited_l_set = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto visited_r_set = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set  = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief get (i,j)-insert graph
         */
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            /**
             * @brief skip edges are not in the current (k-1,k-1)-core
             */
            if((previous_inserted_l_map->count(l) || left_index_map->at(l)->count(k - 1, k - 1))
               && (previous_inserted_r_map->count(r) || right_index_map->at(r)->count(k - 1, k - 1))){
                /**
                 * @brief l is not in the previous (i,j)-core
                 */
                if(!left_index_map->at(l)->count(k, k)){
                    if(!l_map->count(l)){
                        auto core_degree = previous_inserted_l_map->count(l) ? previous_inserted_l_map->at(l):  compute_left_vertex_core_degree(B, right_index_map, previous_inserted_r_map, l, k - 1, k - 1);
                        l_map->insert({l, core_degree});
                    }

                    if(l_map->at(l) >= k){
                        if(visited_l_set->count(l)){
                            ++visited_l_set->at(l);
                        }else{
                            visited_l_set->insert({l, 1});
                        }
                    }
                }

                /**
                 * @brief r is not in the previous (i,j)-core
                 */
                if(!right_index_map->at(r)->count(k, k)){
                    if(!r_map->count(r)){
                        auto core_degree = previous_inserted_r_map->count(r) ? previous_inserted_r_map->at(r): compute_right_vertex_core_degree(B, left_index_map, previous_inserted_l_map, r, k - 1, k - 1);
                        r_map->insert({r, core_degree});
                    }

                    if(r_map->at(r) >= k){
                        if(visited_r_set->count(r)){
                            ++visited_r_set->at(r);
                        }else{
                            visited_r_set->insert({r, 1});
                        }
                    }
                }
            }
        }

        while(!visited_l_set->empty() || !visited_r_set->empty()){
            while(!visited_l_set->empty()){
                auto [l, l_count] = *visited_l_set->begin();
                visited_l_set->erase(l);

                if(l_count == 0 || evicted_l_set->count(l)){
                    continue;
                }

                uint32_t degree = 0;
                for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(!(right_index_map->at(r)->count(k - 1, k - 1) || previous_inserted_r_map->count(r))){
                        continue;
                    }
                    if(inserted_r_map->count(r) || right_index_map->at(r)->count(k, k)){
                        ++degree;
                    }else if(!evicted_r_set->count(r)){
                        if (!r_map->count(r)) {
                            auto core_degree = previous_inserted_r_map->count(r) ?
                                               previous_inserted_r_map->at(r) : compute_right_vertex_core_degree(
                                            B, left_index_map, previous_inserted_l_map, r, k - 1,
                                            k - 1);
                            r_map->insert({r, core_degree});
                        }

                        if (r_map->at(r) >= k) {
                            ++degree;
                            if(visited_r_set->count(r)){
                                ++visited_r_set->at(r);
                            }else{
                                visited_r_set->insert({r, 1});
                            }
                        }
                    }
                }

                if(degree >= k){
                    inserted_l_map->insert({l,degree});
                } else
                {
                    remove_left_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                       visited_l_set, visited_r_set, l, k,k);
                }
            }



            while(!visited_r_set->empty()){
                auto [r, r_count] = *visited_r_set->begin();
                visited_r_set->erase(r);

                if(r_count == 0 || evicted_r_set->count(r)){
                    continue;
                }

                uint32_t degree = 0;
                for(const auto& [l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(!(left_index_map->at(l)->count(k - 1, k - 1) || previous_inserted_l_map->count(l))){
                        continue;
                    }
                    if(inserted_l_map->count(l) || left_index_map->at(l)->count(k, k)){
                        ++degree;
                    }else if(!evicted_l_set->count(l)){
                        if (!l_map->count(l)) {
                            auto core_degree = previous_inserted_l_map->count(l) ?
                                               previous_inserted_l_map->at(l) : compute_left_vertex_core_degree(
                                            B, right_index_map, previous_inserted_r_map, l, k - 1,
                                            k - 1);
                            l_map->insert({l, core_degree});
                        }

                        if (l_map->at(l) >= k) {
                            ++degree;
                            if(visited_l_set->count(l)){
                                ++visited_l_set->at(l);
                            }else{
                                visited_l_set->insert({l, 1});
                            }
                        } else {
                            evicted_l_set->insert(l);
                        }
                    }
                }
                if (degree >= k) {
                    inserted_r_map->insert({r, degree});
                } else {
                    remove_right_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                        visited_l_set, visited_r_set, r, k, k);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::middle_candidate_graph(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                                   const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_order_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_rem_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                   uint32_t k){

        if (!k_order_map->count(k - 1)) {
            k_order_map->insert({k - 1, make_shared<extend_list<int, uint32_t>>()});
        }

        auto order_list = k_order_map->at(k - 1);

//        for (auto p = order_list->get_head(); p; p = p->get_next()) {
//            printf("%u ", p->get_value());
//        }
//        printf("\n");

        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if((!left_index_map->at(l)->count(k - 1, k - 1) && !inserted_l_set->count(l)) ||
               (!right_index_map->at(r)->count(k - 1, k - 1) && !inserted_r_set->count(r)))
            {
                continue;
            }

            if (left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)
                && !(left_index_map->at(l)->count(k + 1, k + 1) && right_index_map->at(r)->count(k + 1, k + 1))) {

                if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k + 1, k + 1)) {
                    ++k_core_degree_map->at(l);
                }

                if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k + 1, k + 1)) {
                    ++k_core_degree_map->at(r);
                }
                continue;
            }

            if (order_list->count_value(l)) {
                auto key = order_list->find_key(l).value();
                for (const auto &[r1, e1]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (order_list->count_value(r1) && order_list->find_key(r1).value() < key) {
                        --k_rem_map->at(r1);
                    }
                }

                order_list->remove(l);
                inserted_l_set->insert(l);
            }

            if(order_list->count_value(r)){
                auto key = order_list->find_key(r).value();
                for (const auto &[l1, e1]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (order_list->count_value(l1) && order_list->find_key(l1).value() < key) {
                        --k_rem_map->at(l1);
                    }
                }

                order_list->remove(r);
                inserted_r_set->insert(r);
            }
        }

        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();


        for (const auto &l: *inserted_l_set) {
            vertex_degree_map->insert({l, k_core_degree_map->at(l)});
            if (k_core_degree_map->at(l) == k - 1) {
                evicted_l_set->insert(l);
            }
        }

        for (const auto &r: *inserted_r_set) {
            vertex_degree_map->insert({r, k_core_degree_map->at(r)});
            if (k_core_degree_map->at(r) == k - 1) {
                evicted_r_set->insert(r);
            }
        }

        auto insertion_vertex_vector = make_shared<vector<uint32_t>>();

        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                insertion_vertex_vector->push_back(l);
                inserted_l_set->erase(l);
                vertex_degree_map->erase(l);

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (vertex_degree_map->count(r)) {
                        --vertex_degree_map->at(r);
                        if (vertex_degree_map->at(r) == k - 1) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                insertion_vertex_vector->push_back(r);
                inserted_r_set->erase(r);
                vertex_degree_map->erase(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (vertex_degree_map->count(l)) {
                        --vertex_degree_map->at(l);
                        if (vertex_degree_map->at(l) == k - 1) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }


//        for (auto p = order_list->get_head(); p; p = p->get_next()) {
//            printf("%u ", p->get_value());
//        }
//        printf("\n");

        auto ext = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto affected_vertex_map = make_shared<map<int, uint32_t>>();

        for (const auto &l: *inserted_l_set) {
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (order_list->count_value(r)) {
                    if (!ext->count(r)) {
                        ext->insert({r, 0});
                        affected_vertex_map->insert({order_list->find_key(r).value(), r});
                    }
                    ++ext->at(r);
                }
            }
        }

        for (const auto &r: *inserted_r_set) {
            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if (order_list->count_value(l)) {
                    if (!ext->count(l)) {
                        ext->insert({l, 0});
                        affected_vertex_map->insert({order_list->find_key(l).value(), l});
                    }
                    ++ext->at(l);
                }
            }
        }

        for(auto iter = insertion_vertex_vector->rbegin(); iter!=insertion_vertex_vector->rend(); ++iter){
            order_list->left_insert(*iter);
        }

        while (!affected_vertex_map->empty()) {
            auto [key, v] = *affected_vertex_map->begin();
            affected_vertex_map->erase(key);

            if (B->get_left_vertex(v)) {
                auto &l = v;
                if (ext->at(l) + k_rem_map->at(l) >= k) {
                    order_list->remove(l);
                    inserted_l_set->insert(l);
                    vertex_degree_map->insert({l, ext->at(l) + k_rem_map->at(l)});
                    ext->at(l) = 0;

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (order_list->count_value(r) && order_list->find_key(r).value() > key) {
                            if (!ext->count(r)) {
                                ext->insert({r, 0});
                                affected_vertex_map->insert({order_list->find_key(r).value(), r});
                            }
                            ++ext->at(r);
                        }
                    }
                }else {
                    k_rem_map->at(l) = ext->at(l) + k_rem_map->at(l);
                    ext->at(l) = 0;

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (vertex_degree_map->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) == k - 1) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }

                    remove_vertex(B, l, evicted_l_set, evicted_r_set, inserted_l_set, inserted_r_set, vertex_degree_map,
                                  order_list, k_rem_map,  ext, affected_vertex_map, k, k);

                }
            } else {
                auto &r = v;
                if (ext->at(r) + k_rem_map->at(r) >= k) {
                    order_list->remove(r);
                    inserted_r_set->insert(r);
                    vertex_degree_map->insert({r, ext->at(r) + k_rem_map->at(r)});
                    ext->at(r) = 0;

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (order_list->count_value(l) && order_list->find_key(l).value() > key) {
                            if (!ext->count(l)) {
                                ext->insert({l, 0});
                                affected_vertex_map->insert({order_list->find_key(l).value(), l});
                            }
                            ++ext->at(l);
                        }
                    }
                } else {
                    k_rem_map->at(r) = ext->at(r) + k_rem_map->at(r);
                    ext->at(r) = 0;

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (vertex_degree_map->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) == k - 1) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                    remove_vertex(B, r, evicted_l_set, evicted_r_set, inserted_l_set, inserted_r_set, vertex_degree_map,
                                  order_list, k_rem_map, ext, affected_vertex_map, k, k);
                }
            }
        }
        order_list->reset_order();

        /**
         * @brief update the core degree of updated vertices
         */
        for (const auto &l: *inserted_l_set) {
            k_core_degree_map->at(l) = vertex_degree_map->at(l);
        }

        for (const auto &r: *inserted_r_set) {
            k_core_degree_map->at(r) = vertex_degree_map->at(r);
        }
    }

    void branch_bipartite_core_maintenance::middle_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
            uint32_t k,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
            const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map){
        auto l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto visited_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto visited_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set  = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief get (i,j)-insert graph
         */
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            /**
             * @brief skip edges are not in the current (k-1,k-1)-core
             */
            if((previous_inserted_l_map->count(l) || left_index_map->at(l)->count(k - 1, k - 1))
               && (previous_inserted_r_map->count(r) || right_index_map->at(r)->count(k - 1, k - 1))){
                /**
                 * @brief l is not in the previous (i,j)-core
                 */
                if(!left_index_map->at(l)->count(k, k)){
                    if(!l_map->count(l)){
                        auto core_degree = previous_inserted_l_map->count(l) ?
                                           previous_inserted_l_map->at(l):  compute_left_vertex_core_degree(left_core_degree_map, right_index_map, previous_inserted_r_map, l, k - 1, k - 1);
                        l_map->insert({l, core_degree});
                    }

                    if(l_map->at(l) >= k){
                        if(visited_l_map->count(l)){
                            ++visited_l_map->at(l);
                        }else{
                            visited_l_map->insert({l, 1});
                        }
                    }else{
                        evicted_l_set->insert(l);
                    }
                }

                /**
                 * @brief r is not in the previous (i,j)-core
                 */
                if(!right_index_map->at(r)->count(k, k)){
                    if(!r_map->count(r)){
                        auto core_degree = previous_inserted_r_map->count(r) ?
                                           previous_inserted_r_map->at(r): compute_right_vertex_core_degree(right_core_degree_map, left_index_map, previous_inserted_l_map, r, k - 1, k - 1);
                        r_map->insert({r, core_degree});
                    }

                    if(r_map->at(r) >= k){
                        if(visited_r_map->count(r)){
                            ++visited_r_map->at(r);
                        }else{
                            visited_r_map->insert({r, 1});
                        }
                    }else{
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!visited_l_map->empty() || !visited_r_map->empty()){
            while(!visited_l_map->empty()){
                auto [l, l_count] = *visited_l_map->begin();
                visited_l_map->erase(l);

                if(l_count == 0 || evicted_l_set->count(l)){
                    continue;
                }

                uint32_t degree = 0;
                auto degree_map = left_core_degree_map->at(l);
                for(auto iter = degree_map->lower_bound(k); iter!=degree_map->end();++iter){
                    for(const auto&r:*iter->second){
                        if(!(right_index_map->at(r)->count(k - 1, k - 1) || previous_inserted_r_map->count(r))){
                            continue;
                        }
                        if(inserted_r_map->count(r) || right_index_map->at(r)->count(k, k)){
                            ++degree;
                        }else if(!evicted_r_set->count(r)){
                            if (!r_map->count(r)) {
                                auto core_degree = previous_inserted_r_map->count(r) ?
                                                   previous_inserted_r_map->at(r) : compute_right_vertex_core_degree(
                                                right_core_degree_map, left_index_map, previous_inserted_l_map, r, k - 1,
                                                k - 1);
                                r_map->insert({r, core_degree});
                            }

                            if (r_map->at(r) >= k) {
                                ++degree;
                                if(visited_r_map->count(r)){
                                    ++visited_r_map->at(r);
                                }else{
                                    visited_r_map->insert({r, 1});
                                }
                            }
                        }
                    }
                }

                if(degree >= k){
                    inserted_l_map->insert({l,degree});
                } else
                {
                    remove_left_vertex(left_core_degree_map, right_core_degree_map,
                                       inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                       visited_l_map, visited_r_map,
                                       l, k, k);
                }
            }

            while(!visited_r_map->empty()){
                auto [r, r_count] = *visited_r_map->begin();
                visited_r_map->erase(r);

                if(r_count == 0 || evicted_r_set->count(r)){
                    continue;
                }

                uint32_t degree = 0;
                auto degree_map = right_core_degree_map->at(r);
                for(auto iter = degree_map->lower_bound(k); iter !=degree_map->end();++iter){
                    for(const auto &l:*iter->second){
                        if(!(left_index_map->at(l)->count(k - 1, k - 1) || previous_inserted_l_map->count(l))){
                            continue;
                        }
                        if(inserted_l_map->count(l) || left_index_map->at(l)->count(k, k)){
                            ++degree;
                        }else if(!evicted_l_set->count(l)){
                            if (!l_map->count(l)) {
                                auto core_degree = previous_inserted_l_map->count(l) ?
                                                   previous_inserted_l_map->at(l) : compute_left_vertex_core_degree(
                                                left_core_degree_map, right_index_map, previous_inserted_r_map, l, k - 1,
                                                k - 1);
                                l_map->insert({l, core_degree});
                            }

                            if (l_map->at(l) >= k) {
                                ++degree;
                                if(visited_l_map->count(l)){
                                    ++visited_l_map->at(l);
                                }else{
                                    visited_l_map->insert({l, 1});
                                }
                            }else{
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
                if(degree >= k){
                    inserted_r_map->insert({r,degree});
                } else
                {
                    remove_right_vertex(left_core_degree_map, right_core_degree_map, inserted_l_map, inserted_r_map,
                                        evicted_l_set, evicted_r_set, visited_l_map, visited_r_map,
                                        r, k, k);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::right_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                                                                  uint32_t k,
                                                                  uint32_t j,
                                                                  const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
                                                                  const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map){
        auto r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto visited_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto visited_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set  = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief get (i,j)-insert graph
         */
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            /**
             * @brief skip edges are not in the current (k,j-1)-core
             */
            if((previous_inserted_l_map->count(l) || left_index_map->at(l)->count(k, j - 1))
               &&(previous_inserted_r_map->count(r) || right_index_map->at(r)->count(j - 1, k)))
            {
                /**
                 * @brief l is not in the previous (k,j)-core
                 */
                if(!left_index_map->at(l)->count(k, j))
                {
                    if(visited_l_map->count(l)){
                        ++visited_l_map->at(l);
                    }else{
                        visited_l_map->insert({l, 1});
                    }
                }

                /**
                 * @brief r is not in the previous (k,j)-core
                 */
                if(!right_index_map->at(r)->count(j, k)){
                    if(!r_map->count(r)){
                        auto core_degree = previous_inserted_r_map->count(r) ? previous_inserted_r_map->at(r) : compute_right_vertex_core_degree(B, left_index_map, previous_inserted_l_map, r, k, j - 1);
                        r_map->insert({r, core_degree});
                    }

                    if(r_map->at(r) >= j){
                        if(visited_r_map->count(r)){
                            ++visited_r_map->at(r);
                        }else{
                            visited_r_map->insert({r, 1});
                        }
                    }
                }
            }
        }


        while(!visited_l_map->empty() || !visited_r_map->empty()){
            while(!visited_l_map->empty()){
                auto [l,l_count] = *visited_l_map->begin();
                visited_l_map->erase(l);

                if(l_count == 0 || evicted_l_set->count(l)){
                    continue;
                }

                uint32_t degree = 0;
                for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(inserted_r_map->count(r) || right_index_map->at(r)->count(j, k)){
                        ++degree;
                    }else if((right_index_map->at(r)->count(j - 1, k) || previous_inserted_r_map->count(r)) && !evicted_r_set->count(r)) {
                        if (!r_map->count(r)) {
                            auto core_degree = previous_inserted_r_map->count(r) ?
                                               previous_inserted_r_map->at(r): compute_right_vertex_core_degree(B,left_index_map,
                                                                                                                previous_inserted_l_map,
                                                                                                                r,
                                                                                                                k,
                                                                                                                j -1);
                            r_map->insert({r, core_degree});
                        }

                        if (r_map->at(r) >= j) {
                            ++degree;
                            if(visited_r_map->count(r)){
                                ++visited_r_map->at(r);
                            }else{
                                visited_r_map->insert({r, 1});
                            }
                        }
                    }
                }

                if(degree >= k){
                    inserted_l_map->insert({l, degree});
                } else
                {
                    remove_left_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                       visited_l_map, visited_r_map, l, k, j);
                }
            }

            while(!visited_r_map->empty()){
                auto [r, r_count] = *visited_r_map->begin();
                visited_r_map->erase(r);

                if(r_count == 0 || evicted_r_set->count(r)){
                    continue;
                }

                uint32_t degree = 0;
                for(const auto& [l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    auto max_j = left_index_map->at(l)->get_j(k);
                    if(inserted_l_map->count(l) || max_j > j -1) {
                        ++degree;
                    }else if((max_j == j - 1 || previous_inserted_l_map->count(l)) && !evicted_l_set->count(l)){
                        ++degree;
                        if(visited_l_map->count(l)){
                            ++visited_l_map->at(l);
                        }else{
                            visited_l_map->insert({l, 1});
                        }
                    }
                }
                if (degree >= j) {
                    inserted_r_map->insert({r, degree});
                } else {
                    remove_right_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                        visited_l_map, visited_r_map, r, k, j);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::right_candidate_graph(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                  const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                                  const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                  uint32_t k,
                                                                  uint32_t j) {
        if (!k_core_order_map->count(j - 1)) {
            k_core_order_map->insert({j - 1, make_shared<extend_list<int, uint32_t>>()});
        }
        auto order_list = k_core_order_map->at(j - 1);

//        for(auto p = order_list->get_head();p;p=p->get_next()){
//            printf("%u ", p->get_value());
//        }
//        printf("\n");

        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if ((!left_index_map->at(l)->count(k, j - 1) && !inserted_l_set->count(l)) ||
                (!right_index_map->at(r)->count(j - 1, k) && !inserted_r_set->count(r))) {
                continue;
            }

            if (left_index_map->at(l)->count(k, j) && right_index_map->at(r)->count(j, k)
                && !(left_index_map->at(l)->count(k, j + 1) && right_index_map->at(r)->count(j + 1, k))) {

                if (left_index_map->at(l)->count(k, j) && !(left_index_map->at(l)->count(k, j + 1))) {
                    ++k_core_degree_map->at(l);
                }

                if (right_index_map->at(r)->count(j, k) && !(right_index_map->at(r)->count(j + 1, k))) {
                    ++k_core_degree_map->at(r);
                }
            }

            if (order_list->count_value(l)) {
                auto key = order_list->find_key(l).value();
                for (const auto &[r1, e1]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (order_list->count_value(r1) && order_list->find_key(r1).value() < key) {
                        --k_core_rem_degree_map->at(r1);
                    }
                }

                order_list->remove(l);
                inserted_l_set->insert(l);
            }

            if(order_list->count_value(r)){
                auto key = order_list->find_key(r).value();
                for (const auto &[l1, e1]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (order_list->count_value(l1) && order_list->find_key(l1).value() < key) {
                        --k_core_rem_degree_map->at(l1);
                    }
                }

                order_list->remove(r);
                inserted_r_set->insert(r);
            }
        }

        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &l: *inserted_l_set) {
            vertex_degree_map->insert({l, k_core_degree_map->at(l)});
        }

        for (const auto &r: *inserted_r_set) {
            vertex_degree_map->insert({r, k_core_degree_map->at(r)});
            if (k_core_degree_map->at(r) == j - 1) {
                evicted_r_set->insert(r);
            }
        }

        auto insertion_vertex_vector = make_shared<vector<uint32_t>>();

        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                insertion_vertex_vector->push_back(l);
                inserted_l_set->erase(l);
                vertex_degree_map->erase(l);

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (vertex_degree_map->count(r)) {
                        --vertex_degree_map->at(r);
                        if (vertex_degree_map->at(r) == j - 1) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                insertion_vertex_vector->push_back(r);
                inserted_r_set->erase(r);
                vertex_degree_map->erase(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (vertex_degree_map->count(l)) {
                        --vertex_degree_map->at(l);
                        if (vertex_degree_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }

//        for(auto p = order_list->get_head();p;p=p->get_next()){
//            printf("%u ", p->get_value());
//        }
//        printf("\n");

        auto ext = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto affected_vertex_map = make_shared<map<int, uint32_t>>();
        for (const auto &l: *inserted_l_set) {
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (order_list->count_value(r)) {
                    if (!ext->count(r)) {
                        ext->insert({r, 0});
                        affected_vertex_map->insert({order_list->find_key(r).value(), r});
                    }
                    ++ext->at(r);
                }
            }
        }

        for (const auto &r: *inserted_r_set) {
            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if (order_list->count_value(l)) {
                    if (!ext->count(l)) {
                        ext->insert({l, 0});
                        affected_vertex_map->insert({order_list->find_key(l).value(), l});
                    }
                    ++ext->at(l);
                }
            }
        }

        for(auto iter = insertion_vertex_vector->rbegin(); iter!=insertion_vertex_vector->rend(); ++iter){
            order_list->left_insert(*iter);
        }

        while (!affected_vertex_map->empty()) {
            auto [key, v] = *affected_vertex_map->begin();
            affected_vertex_map->erase(key);

            if (B->get_right_vertex(v)) {
                auto &r = v;
                if (ext->at(r) + k_core_rem_degree_map->at(r) >= j) {
                    order_list->remove(r);
                    inserted_r_set->insert(r);
                    vertex_degree_map->insert({r, ext->at(r) + k_core_rem_degree_map->at(r)});
                    ext->at(r) = 0;

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (order_list->count_value(l) && order_list->find_key(l).value() > key) {
                            if (!ext->count(l)) {
                                ext->insert({l, 0});
                                affected_vertex_map->insert({order_list->find_key(l).value(), l});
                            }
                            ++ext->at(l);
                        }
                    }
                } else if (ext->at(r) + k_core_rem_degree_map->at(r) < j) {
                    k_core_rem_degree_map->at(r) = ext->at(r) + k_core_rem_degree_map->at(r);
                    ext->at(r) = 0;

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (vertex_degree_map->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) < k) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }

                    remove_vertex(B, r, evicted_l_set, evicted_r_set, inserted_l_set, inserted_r_set, vertex_degree_map,
                                  order_list, k_core_rem_degree_map, ext, affected_vertex_map, k, j);
                }
            } else {
                auto &l = v;

                if (ext->at(l) + k_core_rem_degree_map->at(l) >= k) {
                    order_list->remove(l);
                    inserted_l_set->insert(l);
                    vertex_degree_map->insert({l, ext->at(l) + k_core_rem_degree_map->at(l)});
                    ext->at(l) = 0;

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (order_list->count_value(r) && order_list->find_key(r).value() > key) {
                            if (!ext->count(r)) {
                                ext->insert({r, 0});
                                affected_vertex_map->insert({order_list->find_key(r).value(), r});
                            }
                            ++ext->at(r);
                        }
                    }
                }else{
                    k_core_rem_degree_map->at(l) = ext->at(l) + k_core_rem_degree_map->at(l);
                    ext->at(l) = 0;

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (vertex_degree_map->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) == j - 1) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }

                    remove_vertex(B, l, evicted_l_set, evicted_r_set, inserted_l_set, inserted_r_set, vertex_degree_map,
                                  order_list, k_core_rem_degree_map, ext, affected_vertex_map, k, j);
                }
            }
        }

        order_list->reset_order();

        /**
         * @brief update the core degree of updated vertices
         */
        for (const auto &l: *inserted_l_set) {
            k_core_degree_map->at(l) = vertex_degree_map->at(l);
        }

        for (const auto &r: *inserted_r_set) {
            k_core_degree_map->at(r) = vertex_degree_map->at(r);
        }
    }

    void branch_bipartite_core_maintenance::right_candidate_graph(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_vertex_index_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
            uint32_t k,
            uint32_t j,
            const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
            const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_r_map){
        auto r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto visited_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto visited_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set  = make_shared<unordered_set<uint32_t>>();
        /**
         * @brief get (i,j)-insert graph
         */
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            /**
             * @brief skip edges are not in the current (k,j-1)-core
             */
            if((previous_inserted_l_map->count(l) || left_index_map->at(l)->count(k, j - 1))
               && (previous_inserted_r_map->count(r) || right_vertex_index_map->at(r)->count(j - 1, k)))
            {
                /**
                * @brief l is not in the previous (k,j)-core
                */
                if(!left_index_map->at(l)->count(k, j))
                {
                    if(visited_l_map->count(l)){
                        ++visited_l_map->at(l);
                    }else{
                        visited_l_map->insert({l, 1});
                    }
                }

                /**
                 * @brief r is not in the previous (k,j)-core
                 */
                if(!right_vertex_index_map->at(r)->count(j, k)){
                    if(!r_map->count(r)){
                        auto core_degree = previous_inserted_r_map->count(r) ?
                                           previous_inserted_r_map->at(r) : compute_right_vertex_core_degree(right_core_degree_map, left_index_map, previous_inserted_l_map, r, k, j - 1);
                        r_map->insert({r, core_degree});
                    }

                    if(r_map->at(r) >= j){
                        if(visited_r_map->count(r)){
                            ++visited_r_map->at(r);
                        }else{
                            visited_r_map->insert({r, 1});
                        }
                    }else{
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!visited_l_map->empty() || !visited_r_map->empty()){
            while(!visited_l_map->empty()){
                auto [l, l_count] = *visited_l_map->begin();
                visited_l_map->erase(l);

                if(l_count == 0 || evicted_l_set->count(l)){
                    continue;
                }

                uint32_t degree = 0;
                auto degree_map = left_core_degree_map->at(l);
                for(auto iter = degree_map->lower_bound(j); iter!=degree_map->end();++iter){
                    for(const auto &r:*iter->second){
                        if(!(right_vertex_index_map->at(r)->count(j - 1, k) || previous_inserted_r_map->count(r))){
                            continue;
                        }
                        if(inserted_r_map->count(r) || right_vertex_index_map->at(r)->count(j, k)){
                            ++degree;
                        }else if(!evicted_r_set->count(r)){
                            if (!r_map->count(r)) {
                                auto core_degree = previous_inserted_r_map->count(r) ?
                                                   previous_inserted_r_map->at(r) : compute_right_vertex_core_degree(
                                                right_core_degree_map, left_index_map, previous_inserted_l_map, r, k, j - 1);
                                r_map->insert({r, core_degree});
                            }

                            if (r_map->at(r) >= j) {
                                ++degree;
                                if(visited_r_map->count(r)){
                                    ++visited_r_map->at(r);
                                }else{
                                    visited_r_map->insert({r, 1});
                                }
                            }else{
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }

                if(degree >= k){
                    inserted_l_map->insert({l, degree});
                } else
                {
                    remove_left_vertex(left_core_degree_map, right_core_degree_map, inserted_l_map, inserted_r_map,
                                       evicted_l_set, evicted_r_set, visited_l_map, visited_r_map, l, k, j);
                }
            }


            while(!visited_r_map->empty()){
                auto [r, r_count] = *visited_r_map->begin();
                visited_r_map->erase(r);

                if(r_count == 0 || evicted_r_set->count(r)){
                    continue;
                }

                uint32_t degree = 0;
                auto degree_map = right_core_degree_map->at(r);
                for(auto iter = degree_map->lower_bound(k); iter != degree_map->end(); ++iter){
                    for(const auto &l:*iter->second){
                        if(!(left_index_map->at(l)->count(k, j - 1) || previous_inserted_l_map->count(l)) ){
                            continue;
                        }
                        if(inserted_l_map->count(l) || left_index_map->at(l)->count(k, j)){
                            ++degree;
                        }else if(!evicted_l_set->count(l)){
                            ++degree;
                            if (visited_l_map->count(l)) {
                                ++visited_l_map->at(l);
                            } else {
                                visited_l_map->insert({l, 1});
                            }
                        }
                    }
                }
                if(degree >= j){
                    inserted_r_map->insert({r, degree});
                } else
                {
                    remove_right_vertex(left_core_degree_map, right_core_degree_map, inserted_l_map, inserted_r_map,
                                        evicted_l_set, evicted_r_set, visited_l_map, visited_r_map, r, k, j);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::left_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                            uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        while (!inserted_l_map->empty()) {
            uint32_t  i = inserted_l_map->begin()->second;

            for (const auto&[l, l_degree]: *inserted_l_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == i){
                    evicted_l_set->insert(l);
                }
            }
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                inserted_l_map->erase(l);

                for (uint32_t index = k; index <= i; ++index) {
                    new_left_index_map->at(l)->insert(index, k);
                }

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (inserted_r_map->count(r) && inserted_r_map->at(r) >= k) {
                        --inserted_r_map->at(r);
                        if (inserted_r_map->at(r) < k) {
                            inserted_r_map->erase(r);

                            new_right_index_map->at(r)->insert(k, i);

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (inserted_l_map->count(l2) && inserted_l_map->at(l2) > i) {
                                    --inserted_l_map->at(l2);
                                    if (inserted_l_map->at(l2) == i) {
                                        evicted_l_set->insert(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::left_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph>& B,
                                                                            const shared_ptr<mutex>& global_left_mutex,
                                                                            const shared_ptr<mutex>& global_right_mutex,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                            uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        while (!inserted_l_map->empty()) {
            uint32_t  i = inserted_l_map->begin()->second;

            for (const auto&[l, l_degree]: *inserted_l_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == i){
                    evicted_l_set->insert(l);
                }
            }
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                inserted_l_map->erase(l);

                global_left_mutex->lock();
                for (uint32_t index = k; index <= i; ++index) {
                    new_left_index_map->at(l)->insert(index, k);
                }
                global_left_mutex->unlock();

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (inserted_r_map->count(r) && inserted_r_map->at(r) >= k) {
                        --inserted_r_map->at(r);
                        if (inserted_r_map->at(r) < k) {
                            inserted_r_map->erase(r);

                            global_right_mutex->lock();
                            new_right_index_map->at(r)->insert(k, i);
                            global_right_mutex->unlock();

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (inserted_l_map->count(l2) && inserted_l_map->at(l2) > i) {
                                    --inserted_l_map->at(l2);
                                    if (inserted_l_map->at(l2) == i) {
                                        evicted_l_set->insert(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::left_partial_core_decomposition(
            const shared_ptr<abstract_bipartite_graph> &B,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
            const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
            uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for (uint32_t  i = UINT32_MAX; !inserted_l_set->empty(); ++i) {
            for (const auto&l: *inserted_l_set) {
                k_core_degree_map->at(l) = vertex_degree_map->at(l);

                if (vertex_degree_map->at(l) < i) {
                    i = vertex_degree_map->at(l);

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(vertex_degree_map->at(l) == i){
                    evicted_l_set->insert(l);
                }
            }

            for(const auto &r:*inserted_r_set){
                k_core_degree_map->at(r) = vertex_degree_map->at(r);
            }

            if(!k_core_order_map->count(i)){
                k_core_order_map->insert({i, make_shared<extend_list<int, uint32_t>>()});
            }
            auto order_list = k_core_order_map->at(i);

            while(!evicted_l_set->empty() || !evicted_r_set->empty()){
                while (!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);

                    left_mutex_map->at(l)->lock();
                    for (uint32_t index = k; index <= i; ++index) {
                        new_left_index_map->at(l)->insert(index, k);
                    }
                    left_mutex_map->at(l)->unlock();

                    order_list->push_back(l);
                    k_core_rem_degree_map->at(l) = vertex_degree_map->at(l);
                    inserted_l_set->erase(l);
                    vertex_degree_map->erase(l);

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (inserted_r_set->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) < k) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }

                while(!evicted_r_set->empty()) {
                    auto r = *evicted_r_set->begin();
                    evicted_r_set->erase(r);

                    right_mutex_map->at(r)->lock();
                    new_right_index_map->at(r)->insert(k, i);
                    right_mutex_map->at(r)->unlock();

                    order_list->push_back(r);
                    k_core_rem_degree_map->at(r) = vertex_degree_map->at(r);
                    inserted_r_set->erase(r);
                    vertex_degree_map->erase(r);

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (inserted_l_set->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) == i) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::left_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                            const shared_ptr<mutex>& global_left_mutex,
                                                                            const shared_ptr<mutex>& global_right_mutex,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                            uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        while (!inserted_l_map->empty()) {
            uint32_t  i = inserted_l_map->begin()->second;

            for (const auto&[l, l_degree]: *inserted_l_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == i){
                    evicted_l_set->insert(l);
                }
            }
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                inserted_l_map->erase(l);

                global_left_mutex->lock();
                for (uint32_t index = k; index <= i; ++index) {
                    new_left_index_map->at(l)->insert(index, k);
                }
                global_left_mutex->unlock();

                auto r_degree_map = left_core_degree_map->at(l);
                for (auto r_iter = r_degree_map->lower_bound(k); r_iter != r_degree_map->end(); ++r_iter) {
                    for(const auto &r:*r_iter->second){
                        if (inserted_r_map->count(r) && inserted_r_map->at(r) >= k) {
                            --inserted_r_map->at(r);
                            if (inserted_r_map->at(r) < k) {
                                inserted_r_map->erase(r);

                                global_right_mutex->lock();
                                new_right_index_map->at(r)->insert(k, i);
                                global_right_mutex->unlock();

                                auto l_degree_map = right_core_degree_map->at(r);
                                for (auto l_iter = l_degree_map->lower_bound(i); l_iter != l_degree_map->end(); ++l_iter) {
                                    for(const auto &l2:*l_iter->second){
                                        if (inserted_l_map->count(l2) && inserted_l_map->at(l2) > i) {
                                            --inserted_l_map->at(l2);
                                            if (inserted_l_map->at(l2) == i) {
                                                evicted_l_set->insert(l2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void branch_bipartite_core_maintenance::middle_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                              uint32_t k) {
        while (!inserted_l_map->empty() && !inserted_r_map->empty()){
            {
                auto sub_inserted_l_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_l_map);
                auto sub_inserted_r_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_r_map);
                left_partial_core_decomposition(B, sub_inserted_l_map, sub_inserted_r_map, new_left_index_map, new_right_index_map, k);
            }

            {
                auto sub_inserted_l_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_l_map);
                auto sub_inserted_r_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_r_map);
                right_partial_core_decomposition(B, sub_inserted_l_map, sub_inserted_r_map, new_left_index_map, new_right_index_map, k);
            }

            {
                auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
                auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto&[l, l_degree]: *inserted_l_map) {
                    if (l_degree <= k) {
                        evicted_l_set->insert(l);
                    }
                }

                for (const auto &[r, r_degree]: *inserted_r_map) {
                    if (r_degree <= k) {
                        evicted_r_set->insert(r);
                    }
                }


                while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                    while(!evicted_l_set->empty()){
                        auto l = *evicted_l_set->begin();
                        evicted_l_set->erase(l);

                        inserted_l_map->erase(l);

                        for (const auto &[r, e]:*B->get_left_vertex(l)->get_edge_map()) {
                            if (inserted_r_map->count(r) && inserted_r_map->at(r) > k) {
                                --inserted_r_map->at(r);
                                if (inserted_r_map->at(r) == k) {
                                    evicted_r_set->insert(r);
                                }
                            }
                        }
                    }

                    while (!evicted_r_set->empty()) {
                        auto r = *evicted_r_set->begin();
                        evicted_r_set->erase(r);

                        inserted_r_map->erase(r);

                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (inserted_l_map->count(l) && inserted_l_map->at(l) > k) {
                                --inserted_l_map->at(l);
                                if (inserted_l_map->at(l) == k) {
                                    evicted_l_set->insert(l);
                                }
                            }
                        }
                    }
                }
            }
            ++k;
        }
    }

    void branch_bipartite_core_maintenance::middle_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph>&B,
                                                                              const shared_ptr<mutex>& global_left_mutex,
                                                                              const shared_ptr<mutex>& global_right_mutex,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                              uint32_t k,
                                                                              const shared_ptr<thread_pool>& pool) {
        while (!inserted_l_map->empty() && !inserted_r_map->empty()){
            {
                auto sub_inserted_l_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_l_map);
                auto sub_inserted_r_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_r_map);
                pool->submit_task([=]{
                    left_partial_core_decomposition(B, global_left_mutex, global_right_mutex, sub_inserted_l_map, sub_inserted_r_map, new_left_index_map, new_right_index_map, k);
                });

            }

            {
                auto sub_inserted_l_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_l_map);
                auto sub_inserted_r_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_r_map);
                pool->submit_task([=]{
                    right_partial_core_decomposition(B, global_left_mutex, global_right_mutex, sub_inserted_l_map, sub_inserted_r_map, new_left_index_map, new_right_index_map, k);
                });
            }

            {
                auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
                auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto&[l, l_degree]: *inserted_l_map) {
                    if (l_degree <= k) {
                        evicted_l_set->insert(l);
                    }
                }

                for (const auto &[r, r_degree]: *inserted_r_map) {
                    if (r_degree <= k) {
                        evicted_r_set->insert(r);
                    }
                }


                while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                    while(!evicted_l_set->empty()){
                        auto l = *evicted_l_set->begin();
                        evicted_l_set->erase(l);

                        inserted_l_map->erase(l);

                        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (inserted_r_map->count(r) && inserted_r_map->at(r) > k) {
                                --inserted_r_map->at(r);
                                if (inserted_r_map->at(r) == k) {
                                    evicted_r_set->insert(r);
                                }
                            }
                        }
                    }

                    while (!evicted_r_set->empty()) {
                        auto r = *evicted_r_set->begin();
                        evicted_r_set->erase(r);

                        inserted_r_map->erase(r);

                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (inserted_l_map->count(l) && inserted_l_map->at(l) > k) {
                                --inserted_l_map->at(l);
                                if (inserted_l_map->at(l) == k) {
                                    evicted_l_set->insert(l);
                                }
                            }
                        }
                    }
                }
            }
            ++k;
        }
    }

    void branch_bipartite_core_maintenance::middle_partial_core_decomposition(
            const shared_ptr<abstract_bipartite_graph> &B,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
            const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
            const shared_ptr<bipartite_core_order_index> &core_order_index,
            const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
            const shared_ptr<bipartite_core_degree_index> &core_degree_index,
            uint32_t max_k,
            const shared_ptr<thread_pool> &pool) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                inserted_l_set->size() + inserted_r_set->size());
        for (const auto &l: *inserted_l_set) {
            vertex_degree_map->insert({l, core_degree_index->get_middle_map()->at(l)});
        }

        for (const auto &r: *inserted_r_set) {
            vertex_degree_map->insert({r, core_degree_index->get_middle_map()->at(r)});
        }

        auto k_map = make_shared<map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>();

        for (uint32_t k = max_k; true; ++k) {
            /**
             * @brief middle path
             */
            {
                auto k_core_order_map = core_order_index->get_middle_map();
                auto k_core_rem_map = core_rem_degree_index->get_middle_map();
                auto k_core_degree_map = core_degree_index->get_middle_map();

                auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
                auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto &l: *inserted_l_set) {
                    k_core_degree_map->at(l) = vertex_degree_map->at(l);
                    if (vertex_degree_map->at(l) == k) {
                        evicted_l_set->insert(l);
                    }
                }

                for (const auto &r: *inserted_r_set) {
                    k_core_degree_map->at(r) = vertex_degree_map->at(r);
                    if (vertex_degree_map->at(r) == k) {
                        evicted_r_set->insert(r);
                    }
                }

                if (!k_core_order_map->count(k)) {
                    k_core_order_map->insert({k, make_shared<extend_list<int, uint32_t>>()});
                }
                auto order_list = k_core_order_map->at(k);
                while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                    while (!evicted_l_set->empty()) {
                        auto l = *evicted_l_set->begin();
                        evicted_l_set->erase(l);

                        order_list->push_back(l);
                        k_core_rem_map->at(l) = vertex_degree_map->at(l);
                        inserted_l_set->erase(l);
                        vertex_degree_map->erase(l);

                        for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (inserted_r_set->count(r)) {
                                --vertex_degree_map->at(r);
                                if (vertex_degree_map->at(r) == k) {
                                    evicted_r_set->insert(r);
                                }
                            }
                        }
                    }

                    while (!evicted_r_set->empty()) {
                        auto r = *evicted_r_set->begin();
                        evicted_r_set->erase(r);

                        order_list->push_back(r);
                        k_core_rem_map->at(r) = k_core_degree_map->at(r);
                        inserted_r_set->erase(r);
                        vertex_degree_map->erase(r);

                        for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (inserted_l_set->count(l)) {
                                --vertex_degree_map->at(l);
                                if (vertex_degree_map->at(l) == k) {
                                    evicted_l_set->insert(l);
                                }
                            }
                        }
                    }
                }
            }

            if (!vertex_degree_map->empty()) {
                max_k = k;

                k_map->insert({k + 1, make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map)});

                core_order_index->insert_left(k + 1,
                                              make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_left(k + 1, make_shared<unordered_map<uint32_t, uint32_t>>(
                        vertex_degree_map->size()));
                core_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        vertex_degree_map->size()));

                core_order_index->insert_right(k + 1,
                                               make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_right(k + 1,
                                                    make_shared<unordered_map<uint32_t, uint32_t>>(
                                                            vertex_degree_map->size()));
                core_degree_index->insert_right(k + 1, make_shared<unordered_map<uint32_t, uint32_t>>(
                        vertex_degree_map->size()));

                /**
                 * @brief left path
                 */
                {
                    auto left_k_core_order_map = core_order_index->get_left_map(k);
                    auto left_k_core_rem_degree_map = core_rem_degree_index->get_left_map(k);
                    auto left_k_core_degree_map = core_degree_index->get_left_map(k);

                    auto left_sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                            vertex_degree_map);

                    auto right_k_core_order_map = core_order_index->get_right_map(k);
                    auto right_k_core_rem_degree_map = core_rem_degree_index->get_right_map(k);
                    auto right_k_core_degree_map = core_degree_index->get_right_map(k);


                    pool->submit_task([=] {
                        auto right_sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                                left_sub_vertex_degree_map);
                        /**
                        * @brief left path
                        */
                        {
                            auto sub_inserted_l_set = make_shared<unordered_set<uint32_t>>();
                            auto sub_inserted_r_set = make_shared<unordered_set<uint32_t>>();

                            for (const auto &[v, v_degree]: *left_sub_vertex_degree_map) {
                                left_k_core_degree_map->insert({v, v_degree});
                                left_k_core_rem_degree_map->insert({v, 0});
                                if (B->get_left_vertex(v)) {
                                    sub_inserted_l_set->insert(v);
                                } else {
                                    sub_inserted_r_set->insert(v);
                                }
                            }


                            left_partial_core_decomposition(B, left_mutex_map, right_mutex_map,
                                                            sub_inserted_l_set, sub_inserted_r_set,
                                                            left_sub_vertex_degree_map,
                                                            new_left_index_map,
                                                            new_right_index_map,
                                                            left_k_core_order_map,
                                                            left_k_core_rem_degree_map,
                                                            left_k_core_degree_map,
                                                            k);
                        }


                        /**
                         * @brief right path
                         */
                        {
                            auto sub_inserted_l_set = make_shared<unordered_set<uint32_t>>();
                            auto sub_inserted_r_set = make_shared<unordered_set<uint32_t>>();
                            for (const auto &[v, v_degree]: *right_sub_vertex_degree_map) {
                                right_k_core_degree_map->insert({v, 0});
                                right_k_core_rem_degree_map->insert({v, 0});
                                if (B->get_left_vertex(v)) {
                                    sub_inserted_l_set->insert(v);
                                } else {
                                    sub_inserted_r_set->insert(v);
                                }
                            }

                            right_partial_core_decomposition(B, left_mutex_map, right_mutex_map,
                                                             sub_inserted_l_set,
                                                             sub_inserted_r_set,
                                                             right_sub_vertex_degree_map,
                                                             new_left_index_map, new_right_index_map,
                                                             right_k_core_order_map,
                                                             right_k_core_rem_degree_map,
                                                             right_k_core_degree_map,
                                                             k);
                        }


                    });
                }
            } else {
                break;
            }
        }
    }

    void branch_bipartite_core_maintenance::middle_partial_core_decomposition2(
            const shared_ptr<abstract_bipartite_graph> &B,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
            const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
            const shared_ptr<bipartite_core_order_index> &core_order_index,
            const shared_ptr<bipartite_core_rem_degree_index> &core_rem_degree_index,
            const shared_ptr<bipartite_core_degree_index> &core_degree_index,
            uint32_t max_k,
            const shared_ptr<thread_pool> &pool) {
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>(
                inserted_l_set->size() + inserted_r_set->size());
        for (const auto &l: *inserted_l_set) {
            vertex_degree_map->insert({l, core_degree_index->get_middle_map()->at(l)});
        }

        for (const auto &r: *inserted_r_set) {
            vertex_degree_map->insert({r, core_degree_index->get_middle_map()->at(r)});
        }

        auto k_map = make_shared<map<uint32_t, shared_ptr<unordered_map<uint32_t, uint32_t>>>>();
        for (uint32_t k = max_k; true; ++k) {
            auto k_core_order_map = core_order_index->get_middle_map();
            auto k_core_rem_map = core_rem_degree_index->get_middle_map();
            auto k_core_degree_map = core_degree_index->get_middle_map();

            auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
            auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

            for (const auto &l: *inserted_l_set) {
                k_core_degree_map->at(l) = vertex_degree_map->at(l);
                if (vertex_degree_map->at(l) == k) {
                    evicted_l_set->insert(l);
                }
            }

            for (const auto &r: *inserted_r_set) {
                k_core_degree_map->at(r) = vertex_degree_map->at(r);
                if (vertex_degree_map->at(r) == k) {
                    evicted_r_set->insert(r);
                }
            }

            if (!k_core_order_map->count(k)) {
                k_core_order_map->insert({k, make_shared<extend_list<int, uint32_t>>()});
            }
            auto order_list = k_core_order_map->at(k);
            while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                while (!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);

                    order_list->push_back(l);
                    k_core_rem_map->at(l) = vertex_degree_map->at(l);
                    inserted_l_set->erase(l);
                    vertex_degree_map->erase(l);

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (inserted_r_set->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) == k) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }

                while (!evicted_r_set->empty()) {
                    auto r = *evicted_r_set->begin();
                    evicted_r_set->erase(r);

                    order_list->push_back(r);
                    k_core_rem_map->at(r) = k_core_degree_map->at(r);
                    inserted_r_set->erase(r);
                    vertex_degree_map->erase(r);

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (inserted_l_set->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) == k) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }

            if (!vertex_degree_map->empty()) {
                max_k = k;

                k_map->insert({k + 1, make_shared<unordered_map<uint32_t, uint32_t>>(*vertex_degree_map)});

                core_order_index->insert_left(k + 1,
                                              make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_left(k + 1, make_shared<unordered_map<uint32_t, uint32_t>>(
                        vertex_degree_map->size()));
                core_degree_index->insert_left(k, make_shared<unordered_map<uint32_t, uint32_t>>(
                        vertex_degree_map->size()));

                core_order_index->insert_right(k + 1,
                                               make_shared<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>>());
                core_rem_degree_index->insert_right(k + 1,
                                                    make_shared<unordered_map<uint32_t, uint32_t>>(
                                                            vertex_degree_map->size()));
                core_degree_index->insert_right(k + 1, make_shared<unordered_map<uint32_t, uint32_t>>(
                        vertex_degree_map->size()));
            } else {
                break;
            }
        }


        for (const auto &[k, k_vertex_degree_map]: *k_map) {
            pool->submit_task([=] {
                /**
                 * @brief left path
                 */
                {
                    auto k_core_order_map = core_order_index->get_left_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_left_map(k);
                    auto k_core_degree_map = core_degree_index->get_left_map(k);

                    auto sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                            k_vertex_degree_map);

                    auto sub_inserted_l_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_inserted_r_set = make_shared<unordered_set<uint32_t>>();

                    for (const auto &[v, v_degree]: *sub_vertex_degree_map) {
                        k_core_degree_map->insert({v, v_degree});
                        k_core_rem_degree_map->insert({v, 0});
                        if (B->get_left_vertex(v)) {
                            sub_inserted_l_set->insert(v);
                        } else {
                            sub_inserted_r_set->insert(v);
                        }
                    }

                    left_partial_core_decomposition(B, left_mutex_map, right_mutex_map,
                                                    sub_inserted_l_set, sub_inserted_r_set,
                                                    sub_vertex_degree_map,
                                                    new_left_index_map,
                                                    new_right_index_map,
                                                    k_core_order_map,
                                                    k_core_rem_degree_map,
                                                    k_core_degree_map,
                                                    k);
                }
            });
        }

        for (const auto &[k, k_vertex_degree_map]: *k_map) {
            pool->submit_task([=] {
                /**
                 * @brief right path
                 */
                {
                    auto k_core_order_map = core_order_index->get_right_map(k);
                    auto k_core_rem_degree_map = core_rem_degree_index->get_right_map(k);
                    auto k_core_degree_map = core_degree_index->get_right_map(k);

                    auto sub_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                            k_vertex_degree_map);

                    auto sub_inserted_l_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_inserted_r_set = make_shared<unordered_set<uint32_t>>();
                    for (const auto &[v, v_degree]: *sub_vertex_degree_map) {
                        k_core_degree_map->insert({v, 0});
                        k_core_rem_degree_map->insert({v, 0});
                        if (B->get_left_vertex(v)) {
                            sub_inserted_l_set->insert(v);
                        } else {
                            sub_inserted_r_set->insert(v);
                        }
                    }

                    right_partial_core_decomposition(B, left_mutex_map, right_mutex_map,
                                                     sub_inserted_l_set,
                                                     sub_inserted_r_set,
                                                     sub_vertex_degree_map,
                                                     new_left_index_map, new_right_index_map,
                                                     k_core_order_map,
                                                     k_core_rem_degree_map,
                                                     k_core_degree_map,
                                                     k);
                }
            });
        }
    }

    void branch_bipartite_core_maintenance::middle_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                              const shared_ptr<mutex>& global_left_mutex,
                                                                              const shared_ptr<mutex>& global_right_mutex,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                              uint32_t k,
                                                                              const shared_ptr<thread_pool>& pool) {
        while (!inserted_l_map->empty() && !inserted_r_map->empty()){
            {
                auto sub_inserted_l_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_l_map);
                auto sub_inserted_r_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_r_map);
                pool->submit_task([=] {
                    left_partial_core_decomposition(left_core_degree_map, right_core_degree_map, global_left_mutex, global_right_mutex,
                                                    sub_inserted_l_map, sub_inserted_r_map, new_left_index_map,
                                                    new_right_index_map, k);
                });
            }

            {
                auto sub_inserted_l_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_l_map);
                auto sub_inserted_r_map = container_copy::to_unordered_map<uint32_t,uint32_t>(inserted_r_map);
                pool->submit_task([=]{
                    right_partial_core_decomposition(left_core_degree_map, right_core_degree_map, global_left_mutex, global_right_mutex, sub_inserted_l_map, sub_inserted_r_map, new_left_index_map, new_right_index_map, k);
                });
            }

            {
                auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
                auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto&[l, l_degree]: *inserted_l_map) {
                    if (l_degree == k) {
                        evicted_l_set->insert(l);
                    }
                }

                for (const auto &[r, r_degree]: *inserted_r_map) {
                    if (r_degree == k) {
                        evicted_r_set->insert(r);
                    }
                }


                while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                    while(!evicted_l_set->empty()){
                        auto l = *evicted_l_set->begin();
                        evicted_l_set->erase(l);

                        inserted_l_map->erase(l);

                        auto degree_map = left_core_degree_map->at(l);
                        for (auto iter = degree_map->lower_bound(k); iter!=degree_map->end(); ++iter) {
                            for(const auto &r:*iter->second){
                                if (inserted_r_map->count(r) && inserted_r_map->at(r) > k) {
                                    --inserted_r_map->at(r);
                                    if (inserted_r_map->at(r) == k) {
                                        evicted_r_set->insert(r);
                                    }
                                }
                            }
                        }
                    }

                    while (!evicted_r_set->empty()) {
                        auto r = *evicted_r_set->begin();
                        evicted_r_set->erase(r);

                        inserted_r_map->erase(r);

                        auto degree_map = right_core_degree_map->at(r);
                        for (auto iter = degree_map->lower_bound(k); iter!=degree_map->end(); ++iter) {
                            for(const auto&l:*iter->second){
                                if (inserted_l_map->count(l) && inserted_l_map->at(l) > k) {
                                    --inserted_l_map->at(l);
                                    if (inserted_l_map->at(l) == k) {
                                        evicted_l_set->insert(l);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ++k;
        }
    }

    void branch_bipartite_core_maintenance::right_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph> &B,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                             uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        while (!inserted_r_map->empty()) {
            uint32_t j = inserted_r_map->begin()->second;

            for (const auto &[r, r_degree]: *inserted_r_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (r_degree == j) {
                    evicted_r_set->insert(r);
                }
            }
            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                inserted_r_map->erase(r);

                for (uint32_t index = k; index <= j; ++index) {
                    new_right_index_map->at(r)->insert(index, k);
                }

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (inserted_l_map->count(l) && inserted_l_map->at(l) >= k) {
                        --inserted_l_map->at(l);
                        if (inserted_l_map->at(l) < k) {
                            new_left_index_map->at(l)->insert(k, j);

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (inserted_r_map->count(r2) && inserted_r_map->at(r2) > j) {
                                    --inserted_r_map->at(r2);
                                    if (inserted_r_map->at(r2) == j) {
                                        evicted_r_set->insert(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::right_partial_core_decomposition(const shared_ptr<abstract_bipartite_graph>& B,
                                                                             const shared_ptr<mutex>& global_left_mutex,
                                                                             const shared_ptr<mutex>& global_right_mutex,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                             uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        while (!inserted_r_map->empty()) {
            uint32_t j = inserted_r_map->begin()->second;

            for (const auto &[r, r_degree]: *inserted_r_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (r_degree == j) {
                    evicted_r_set->insert(r);
                }
            }
            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                inserted_r_map->erase(r);

                global_right_mutex->lock();
                for (uint32_t index = k; index <= j; ++index) {
                    new_right_index_map->at(r)->insert(index, k);
                }
                global_right_mutex->unlock();

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (inserted_l_map->count(l) && inserted_l_map->at(l) >= k) {
                        --inserted_l_map->at(l);
                        if (inserted_l_map->at(l) < k) {

                            global_left_mutex->lock();
                            new_left_index_map->at(l)->insert(k, j);
                            global_left_mutex->unlock();

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (inserted_r_map->count(r2) && inserted_r_map->at(r2) > j) {
                                    --inserted_r_map->at(r2);
                                    if (inserted_r_map->at(r2) == j) {
                                        evicted_r_set->insert(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::right_partial_core_decomposition(
            const shared_ptr<abstract_bipartite_graph> &B,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
            const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
            const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
            uint32_t k) {

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for (uint32_t j = UINT32_MAX; !inserted_r_set->empty(); ++j) {

            for(const auto &l:*inserted_l_set){
                k_core_degree_map->at(l) = vertex_degree_map->at(l);
            }

            for (const auto &r: *inserted_r_set) {
                k_core_degree_map->at(r) = vertex_degree_map->at(r);

                if (vertex_degree_map->at(r) < j) {
                    j = vertex_degree_map->at(r);

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (vertex_degree_map->at(r) == j) {
                    evicted_r_set->insert(r);
                }
            }

            if(!k_core_order_map->count(j)){
                k_core_order_map->insert({j, make_shared<extend_list<int, uint32_t>>()});
            }
            auto order_list = k_core_order_map->at(j);

            while(!evicted_l_set->empty() || !evicted_r_set->empty()){


                while (!evicted_r_set->empty()) {
                    auto r = *evicted_r_set->begin();
                    evicted_r_set->erase(r);

                    right_mutex_map->at(r)->lock();
                    for (uint32_t index = k; index <= j; ++index) {
                        new_right_index_map->at(r)->insert(index, k);
                    }
                    right_mutex_map->at(r)->unlock();

                    order_list->push_back(r);
                    k_core_rem_degree_map->at(r) = vertex_degree_map->at(r);
                    inserted_r_set->erase(r);

                    for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                        if (inserted_l_set->count(l)) {
                            --vertex_degree_map->at(l);
                            if (vertex_degree_map->at(l) < k) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }

                while(!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);

                    left_mutex_map->at(l)->lock();
                    new_left_index_map->at(l)->insert(k, j);
                    left_mutex_map->at(l)->unlock();

                    order_list->push_back(l);
                    k_core_rem_degree_map->at(l) = vertex_degree_map->at(l);
                    inserted_l_set->erase(l);

                    for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                        if (inserted_r_set->count(r)) {
                            --vertex_degree_map->at(r);
                            if (vertex_degree_map->at(r) == j) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::right_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                             const shared_ptr<mutex>& global_left_mutex,
                                                                             const shared_ptr<mutex>& global_right_mutex,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                             uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        while (!inserted_r_map->empty()) {
            uint32_t j = inserted_r_map->begin()->second;

            for (const auto &[r, r_degree]: *inserted_r_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if (r_degree == j) {
                    evicted_r_set->insert(r);
                }
            }
            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                inserted_r_map->erase(r);

                global_right_mutex->lock();
                for (uint32_t index = k; index <= j; ++index) {
                    new_right_index_map->at(r)->insert(index, k);
                }
                global_right_mutex->unlock();

                auto l_degree_map = right_core_degree_map->at(r);
                for (auto l_iter = l_degree_map->lower_bound(k); l_iter != l_degree_map->end(); ++l_iter) {
                    for(const auto &l:*l_iter->second){
                        if (inserted_l_map->count(l) && inserted_l_map->at(l) >= k) {
                            --inserted_l_map->at(l);
                            if (inserted_l_map->at(l) < k) {

                                global_left_mutex->lock();
                                new_left_index_map->at(l)->insert(k, j);
                                global_left_mutex->unlock();

                                auto r_degree_map = left_core_degree_map->at(l);
                                for (auto r_iter = r_degree_map->lower_bound(j); r_iter != r_degree_map->end(); ++r_iter) {
                                    for(const auto&r2:*r_iter->second){
                                        if (inserted_r_map->count(r2) && inserted_r_map->at(r2) > j) {
                                            --inserted_r_map->at(r2);
                                            if (inserted_r_map->at(r2) == j) {
                                                evicted_r_set->insert(r2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::left_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                      const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                                      const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                                      uint32_t i,
                                                                      uint32_t k,
                                                                      const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                                      const shared_ptr<unordered_set<uint32_t>> &removed_r_set) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        for(const auto&e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(left_index_map->at(l)->count(i, k) && right_index_map->at(r)->count(k, i)){
                if (B->get_left_vertex(l) && (!left_index_map->at(l)->count(i + 1, k) || previous_removed_l_set->count(l))) {
                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, i, k);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) < i) {
                        evicted_l_set->insert(l);
                    }
                }

                if (B->get_right_vertex(r) && (!right_index_map->at(r)->count(k, i + 1) || previous_removed_r_set->count(r))) {
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, i, k);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            for(const auto &l:*evicted_l_set){
                for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(!right_index_map->at(r)->count(k, i) || removed_r_set->count(r)){
                        continue;
                    }

                    if(!right_index_map->at(r)->count(k, i + 1) || previous_removed_r_set->count(r)){
                        if (!r_map->count(r)) {
                            auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, i, k);
                            r_map->insert({r, r_degree});
                        }

                        --r_map->at(r);
                        if (r_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }
            removed_l_set->merge(*evicted_l_set);

            for(const auto r:*evicted_r_set){
                for(const auto&[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(!left_index_map->at(l)->count(i, k) || removed_l_set->count(l)){
                        continue;
                    }

                    if(!left_index_map->at(l)->count( i + 1, k) || previous_removed_l_set->count(l)){
                        if (!l_map->count(l)) {
                            auto l_degree = compute_left_vertex_core_degree(B, right_index_map,  l, i, k);
                            l_map->insert({l, l_degree});
                        }

                        --l_map->at(l);
                        if (l_map->at(l) < i) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
            removed_r_set->merge(*evicted_r_set);
        }
    }

    void branch_bipartite_core_maintenance::left_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                                                      const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                      const shared_ptr<vector<uint32_t>> &previous_removed_vector,
                                                                      uint32_t i,
                                                                      uint32_t k,
                                                                      const shared_ptr<vector<uint32_t>> &removed_vector) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        if (!k_core_order_map->count(i)) {
            k_core_order_map->insert({i, make_shared<extend_list<int, uint32_t>>()});
        }
        auto order_list = k_core_order_map->at(i);

        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if (left_index_map->at(l)->count(i, k) && right_index_map->at(r)->count(k, i)) {

                if(!order_list->count_value(l) && !order_list->count_value(r)){
                    continue;
                }

                if (order_list->count_value(l) && order_list->count_value(r)) {
                    --k_core_degree_map->at(l);
                    if (k_core_degree_map->at(l) < i) {
                        evicted_l_set->insert(l);
                    }

                    --k_core_degree_map->at(r);
                    if (k_core_degree_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    }

                    if (order_list->find(l)->get_key() < order_list->find(r)->get_key()) {
                        --k_core_rem_degree_map->at(l);
                    } else {
                        --k_core_rem_degree_map->at(r);
                    }

                } else if (order_list->count_value(l)) {
                    --k_core_degree_map->at(l);
                    if (k_core_degree_map->at(l) < i) {
                        evicted_l_set->insert(l);
                    } else {
                        --k_core_rem_degree_map->at(l);
                    }
                } else {
                    --k_core_degree_map->at(r);
                    if (k_core_degree_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    } else {
                        --k_core_rem_degree_map->at(r);
                    }
                }
            }
        }

        for (const auto &v: *previous_removed_vector) {
            if(left_index_map->count(v)){
                auto &l = v;
                k_core_degree_map->at(l) = compute_left_vertex_core_degree(B, right_index_map, l, i, k);
                if (k_core_degree_map->at(l) < i) {
                    evicted_l_set->insert(l);
                }
                order_list->push_back(l);
            }else{
                auto &r = v;
                k_core_degree_map->at(r) = compute_right_vertex_core_degree(B, left_index_map, r, i, k);
                if (k_core_degree_map->at(r) < k) {
                    evicted_r_set->insert(r);
                }
                order_list->push_back(r);
            }
        }


        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                auto key = order_list->find(l)->get_key();
                order_list->remove(l);
                removed_vector->push_back(l);
                k_core_rem_degree_map->at(l) = k_core_degree_map->at(l);

                if(!B->get_left_vertex(l)){
                    continue;
                }

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (order_list->count_value(r)) {
                        --k_core_degree_map->at(r);
                        if (k_core_degree_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        } else if (order_list->find(r)->get_key() < key) {
                            --k_core_rem_degree_map->at(r);
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                auto key = order_list->find(r)->get_key();
                order_list->remove(r);
                removed_vector->push_back(r);
                k_core_rem_degree_map->at(r) = k_core_degree_map->at(r);

                if(!B->get_right_vertex(r)){
                    continue;
                }

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (order_list->count_value(l)) {
                        --k_core_degree_map->at(l);
                        if (k_core_degree_map->at(l) < i) {
                            evicted_l_set->insert(l);
                        } else if (order_list->find(l)->get_key() < key) {
                            --k_core_rem_degree_map->at(l);
                        }
                    }
                }
            }
        }

        order_list->reset_order();
    }

    void branch_bipartite_core_maintenance::left_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                      const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                                      const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                                      uint32_t i,
                                                                      uint32_t k,
                                                                      const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                                      const shared_ptr<unordered_set<uint32_t>> &removed_r_set) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(left_index_map->at(l)->count(i, k) && right_index_map->at(r)->count(k, i)){
                if (B->get_left_vertex(l) && (!left_index_map->at(l)->count(i + 1, k) || previous_removed_l_set->count(l))) {
                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(left_core_degree_map, right_index_map , l, i, k);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) < i) {
                        evicted_l_set->insert(l);
                    }
                }

                if (B->get_right_vertex(r) && (!right_index_map->at(r)->count(k, i + 1) || previous_removed_r_set->count(r))) {
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(right_core_degree_map, left_index_map, r, i, k);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){

            for(const auto &l:*evicted_l_set){
                auto degree_map = left_core_degree_map->at(l);
                for(auto iter = degree_map->lower_bound(k); iter !=degree_map->end(); ++iter){
                    for(const auto&r:*iter->second){
                        if(!right_index_map->at(r)->count(k, i) || removed_r_set->count(r)){
                            continue;
                        }
                        if(!right_index_map->at(r)->count(k, i+1) || previous_removed_r_set->count(r)){
                            if (!r_map->count(r)) {
                                auto r_degree = compute_right_vertex_core_degree(right_core_degree_map,left_index_map, r, i, k);
                                r_map->insert({r, r_degree});
                            }

                            --r_map->at(r);
                            if (r_map->at(r) < k) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }
            }
            removed_l_set->merge(*evicted_l_set);

            for(const auto &r:*evicted_r_set){
                auto degree_map = right_core_degree_map->at(r);
                for(auto iter = degree_map->lower_bound(i); iter!=degree_map->end(); ++iter){
                    for(const auto&l:*iter->second){
                        if(!left_index_map->at(l)->count(i, k) || removed_l_set->count(l)){
                            continue;
                        }
                        if(!left_index_map->at(l)->count(i + 1, k) || previous_removed_l_set->count(l)){
                            if (!l_map->count(l)) {
                                auto l_degree = compute_left_vertex_core_degree(left_core_degree_map, right_index_map, l, i, k);
                                l_map->insert({l, l_degree});
                            }

                            --l_map->at(l);
                            if (l_map->at(l) < i) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
            removed_r_set->merge(*evicted_r_set);
        }
    }

    void branch_bipartite_core_maintenance::middle_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                        const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                                        const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                                        uint32_t k,
                                                                        const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                                        const shared_ptr<unordered_set<uint32_t>> &removed_r_set) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)){
                if (B->get_left_vertex(l) && (!left_index_map->at(l)->count(k + 1, k + 1) || previous_removed_l_set->count(l))) {
                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, k, k);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }
                }

                if (B->get_right_vertex(r) && (!right_index_map->at(r)->count(k + 1, k + 1) || previous_removed_r_set->count(r))) {
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, k, k);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            for(const auto &l:*evicted_l_set){
                for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(!right_index_map->at(r)->count(k, k) || removed_r_set->count(r)){
                        continue;
                    }

                    if(!right_index_map->at(r)->count(k + 1, k + 1) || previous_removed_r_set->count(r)){
                        if (!r_map->count(r)) {
                            auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, k, k);
                            r_map->insert({r, r_degree});
                        }

                        --r_map->at(r);
                        if (r_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }
            removed_l_set->merge(*evicted_l_set);

            for(const auto &r:*evicted_r_set){
                for(const auto&[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(!left_index_map->at(l)->count(k, k) || removed_l_set->count(l) ){
                        continue;
                    }

                    if (!left_index_map->at(l)->count(k + 1, k + 1) || previous_removed_l_set->count(l)){
                        if (!l_map->count(l)) {
                            auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, k, k);
                            l_map->insert({l, l_degree});
                        }

                        --l_map->at(l);
                        if (l_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
            removed_r_set->merge(*evicted_r_set);
        }
    }

    void branch_bipartite_core_maintenance::middle_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                        const shared_ptr<vector<uint32_t>> &previous_removed_vector,
                                                                        uint32_t k,
                                                                        const shared_ptr<vector<uint32_t>> &removed_vector) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        if (!k_core_order_map->count(k)) {
            k_core_order_map->insert({k, make_shared<extend_list<int, uint32_t>>()});
        }
        auto order_list = k_core_order_map->at(k);

        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if (left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)) {

                if(!order_list->count_value(l) && !order_list->count_value(r)){
                    continue;
                }

                if (order_list->count_value(l) && order_list->count_value(r)) {
                    --k_core_degree_map->at(l);
                    if (k_core_degree_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }

                    --k_core_degree_map->at(r);
                    if (k_core_degree_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    }

                    if (order_list->find(l)->get_key() < order_list->find(r)->get_key()) {
                        --k_core_rem_degree_map->at(l);
                    } else {
                        --k_core_rem_degree_map->at(r);
                    }

                } else if (!order_list->count_value(r)) {
                    --k_core_degree_map->at(l);
                    if (k_core_degree_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }else{
                        --k_core_rem_degree_map->at(l);
                    }
                } else {
                    --k_core_degree_map->at(r);
                    if (k_core_degree_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    }else{
                        --k_core_rem_degree_map->at(r);
                    }
                }
            }
        }

        for (const auto &v: *previous_removed_vector) {
            if(left_index_map->count(v)){
                auto &l = v;
                k_core_degree_map->at(l) = compute_left_vertex_core_degree(B, right_index_map, l, k, k);
                if (k_core_degree_map->at(l) < k) {
                    evicted_l_set->insert(l);
                }
                order_list->push_back(l);
            }else{
                auto &r = v;
                k_core_degree_map->at(r) = compute_right_vertex_core_degree(B, left_index_map, r, k, k);
                if (k_core_degree_map->at(r) < k) {
                    evicted_r_set->insert(r);
                }
                order_list->push_back(r);
            }
        }


        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                auto key = order_list->find(l)->get_key();
                order_list->remove(l);
                removed_vector->push_back(l);
                k_core_rem_degree_map->at(l) = k_core_degree_map->at(l);

                if(!B->get_left_vertex(l)){
                    continue;
                }

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (order_list->count_value(r)) {
                        --k_core_degree_map->at(r);
                        if (k_core_degree_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        } else if (order_list->find(r)->get_key() < key) {
                            --k_core_rem_degree_map->at(r);
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                auto key = order_list->find(r)->get_key();
                order_list->remove(r);
                removed_vector->push_back(r);
                k_core_rem_degree_map->at(r) = k_core_degree_map->at(r);

                if(!B->get_right_vertex(r)){
                    continue;
                }

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (order_list->count_value(l)) {
                        --k_core_degree_map->at(l);
                        if (k_core_degree_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        } else if (order_list->find(l)->get_key() < key) {
                            --k_core_rem_degree_map->at(l);
                        }
                    }
                }
            }
        }

        order_list->reset_order();
    }

    void branch_bipartite_core_maintenance::middle_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                        const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                                        const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                                        uint32_t k,
                                                                        const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                                        const shared_ptr<unordered_set<uint32_t>> &removed_r_set) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        for(const auto&e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(left_index_map->at(l)->count(k, k) && right_index_map->at(r)->count(k, k)){
                if (B->get_left_vertex(l) && (!left_index_map->at(l)->count(k + 1, k + 1) || previous_removed_l_set->count(l))) {
                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(left_core_degree_map, right_index_map, l, k, k);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }
                }

                if (B->get_right_vertex(r) && (!right_index_map->at(r)->count(k + 1, k + 1) || previous_removed_r_set->count(r))) {
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(right_core_degree_map, left_index_map, r, k, k);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) < k) {
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            for(const auto &l:*evicted_l_set){
                auto degree_map = left_core_degree_map->at(l);
                for(auto iter = degree_map->lower_bound(k); iter!=degree_map->end(); ++iter){
                    for(const auto&r:*iter->second){
                        if(!right_index_map->at(r)->count(k, k) || removed_r_set->count(r)){
                            continue;
                        }
                        if(!right_index_map->at(r)->count(k + 1, k + 1) || previous_removed_r_set->count(r)){
                            if (!r_map->count(r)) {
                                auto r_degree = compute_right_vertex_core_degree(right_core_degree_map, left_index_map, r, k, k);
                                r_map->insert({r, r_degree});
                            }

                            --r_map->at(r);
                            if (r_map->at(r) < k) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }
            }
            removed_l_set->merge(*evicted_l_set);

            for(const auto &r:*evicted_r_set){
                auto degree_map = right_core_degree_map->at(r);
                for(auto iter = degree_map->lower_bound(k); iter!=degree_map->end();++iter){
                    for(const auto&l:*iter->second){
                        if(!left_index_map->at(l)->count(k, k) || removed_l_set->count(l)){
                            continue;
                        }
                        if(!left_index_map->at(l)->count(k + 1, k + 1) || previous_removed_l_set->count(l)){
                            if (!l_map->count(l)) {
                                auto l_degree = compute_left_vertex_core_degree(left_core_degree_map, right_index_map, l, k, k);
                                l_map->insert({l, l_degree});
                            }

                            --l_map->at(l);
                            if (l_map->at(l) < k) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
            removed_r_set->merge(*evicted_r_set);
        }
    }

    void branch_bipartite_core_maintenance::right_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                                       const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                                       uint32_t k,
                                                                       uint32_t j,
                                                                       const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                                       const shared_ptr<unordered_set<uint32_t>> &removed_r_set) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(left_index_map->at(l)->count(k, j) && right_index_map->at(r)->count(j, k)){
                if(B->get_left_vertex(l) && (!left_index_map->at(l)->count(k, j + 1) || previous_removed_l_set->count(l))) {
                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, k, j);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }
                }

                if(B->get_right_vertex(r) && (!right_index_map->at(r)->count(j + 1, k) || previous_removed_r_set->count(r))) {
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, k, j);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) < j) {
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            for(const auto &l:*evicted_l_set){
                for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(!right_index_map->at(r)->count(j, k) || removed_r_set->count(r)){
                        continue;
                    }
                    if(!right_index_map->at(r)->count(j + 1,k) || previous_removed_r_set->count(r)){
                        if (!r_map->count(r)) {
                            auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, k, j);
                            r_map->insert({r, r_degree});
                        }

                        --r_map->at(r);
                        if (r_map->at(r) < j) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }
            removed_l_set->merge(*evicted_l_set);

            for(const auto &r:*evicted_r_set){
                for(const auto&[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(!left_index_map->at(l)->count(k, j) || removed_l_set->count(l)){
                        continue;
                    }
                    if(!left_index_map->at(l)->count(k, j+1) || previous_removed_l_set->count(l)){
                        if (!l_map->count(l)) {
                            auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, k, j);
                            l_map->insert({l, l_degree});
                        }

                        --l_map->at(l);
                        if (l_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
            removed_r_set->merge(*evicted_r_set);
        }
    }

    void branch_bipartite_core_maintenance::right_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, uint32_t>>>> &k_core_order_map,
                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_rem_degree_map,
                                                                       const shared_ptr<unordered_map<uint32_t, uint32_t>> &k_core_degree_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                       const shared_ptr<vector<uint32_t>> &previous_removed_vector,
                                                                       uint32_t k,
                                                                       uint32_t j,
                                                                       const shared_ptr<vector<uint32_t>> &removed_vector) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        if (!k_core_order_map->count(j)) {
            k_core_order_map->insert({j, make_shared<extend_list<int, uint32_t>>()});
        }
        auto order_list = k_core_order_map->at(j);

        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if (left_index_map->at(l)->count(k, j) && right_index_map->at(r)->count(j, k)) {

                if(!order_list->count_value(l) && !order_list->count_value(r)){
                    continue;
                }

                if (order_list->count_value(l) && order_list->count_value(r)) {
                    --k_core_degree_map->at(l);
                    if (k_core_degree_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }

                    --k_core_degree_map->at(r);
                    if (k_core_degree_map->at(r) < j) {
                        evicted_r_set->insert(r);
                    }

                    if (order_list->find(l)->get_key() < order_list->find(r)->get_key()) {
                        --k_core_rem_degree_map->at(l);
                    } else {
                        --k_core_rem_degree_map->at(r);
                    }

                } else if(order_list->count_value(l)) {
                    --k_core_degree_map->at(l);
                    if (k_core_degree_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }else{
                        --k_core_rem_degree_map->at(l);
                    }
                } else{
                    --k_core_degree_map->at(r);
                    if (k_core_degree_map->at(r) < j) {
                        evicted_r_set->insert(r);
                    }else{
                        --k_core_rem_degree_map->at(r);
                    }
                }
            }
        }

        for (const auto &v: *previous_removed_vector) {
            if(left_index_map->count(v)){
                auto &l = v;
                k_core_degree_map->at(l) = compute_left_vertex_core_degree(B, right_index_map, l, k, j);
                if (k_core_degree_map->at(l) < k) {
                    evicted_l_set->insert(l);
                }
                order_list->push_back(l);
            }else{
                auto &r = v;
                k_core_degree_map->at(r) = compute_right_vertex_core_degree(B, left_index_map, r, k, j);
                if (k_core_degree_map->at(r) < j) {
                    evicted_r_set->insert(r);
                }
                order_list->push_back(r);
            }
        }


        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                auto key = order_list->find(l)->get_key();
                order_list->remove(l);
                removed_vector->push_back(l);
                k_core_rem_degree_map->at(l) = k_core_degree_map->at(l);

                if(!B->get_left_vertex(l)){
                    continue;
                }

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (order_list->count_value(r)) {
                        --k_core_degree_map->at(r);
                        if (k_core_degree_map->at(r) < j) {
                            evicted_r_set->insert(r);
                        } else if (order_list->find(r)->get_key() < key) {
                            --k_core_rem_degree_map->at(r);
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                auto key = order_list->find(r)->get_key();
                order_list->remove(r);
                removed_vector->push_back(r);
                k_core_rem_degree_map->at(r) = k_core_degree_map->at(r);

                if(!B->get_right_vertex(r)){
                    continue;
                }

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (order_list->count_value(l)) {
                        --k_core_degree_map->at(l);
                        if (k_core_degree_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        } else if (order_list->find(l)->get_key() < key) {
                            --k_core_rem_degree_map->at(l);
                        }
                    }
                }
            }
        }

        order_list->reset_order();
    }

    void branch_bipartite_core_maintenance::right_removal_partial_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &previous_removed_l_set,
                                                                       const shared_ptr<unordered_set<uint32_t>> &previous_removed_r_set,
                                                                       uint32_t k,
                                                                       uint32_t j,
                                                                       const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                                       const shared_ptr<unordered_set<uint32_t>> &removed_r_set) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(left_index_map->at(l)->count(k, j) && right_index_map->at(r)->count(j, k)){
                if (B->get_left_vertex(l) && (!left_index_map->at(l)->count(k, j + 1) || previous_removed_l_set->count(l))) {
                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(left_core_degree_map, right_index_map, l, k, j);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) < k) {
                        evicted_l_set->insert(l);
                    }
                }

                if (B->get_right_vertex(r) && (!right_index_map->at(r)->count(j + 1, k) || previous_removed_r_set->count(r))) {
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(right_core_degree_map, left_index_map, r, k, j);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) < j) {
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            for(const auto &l:*evicted_l_set){
                auto degree_map = left_core_degree_map->at(l);
                for(auto iter = degree_map->lower_bound(j); iter!=degree_map->end(); ++iter){
                    for(const auto&r:*iter->second){
                        if(!right_index_map->at(r)->count(j, k) || removed_r_set->count(r)){
                            continue;
                        }
                        if(!right_index_map->at(r)->count(j + 1,k) || previous_removed_r_set->count(r)){
                            if (!r_map->count(r)) {
                                auto r_degree = compute_right_vertex_core_degree(right_core_degree_map, left_index_map, r, k, j);
                                r_map->insert({r, r_degree});
                            }

                            --r_map->at(r);
                            if (r_map->at(r) < j) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }
            }
            removed_l_set->merge(*evicted_l_set);

            for(const auto &r:*evicted_r_set){
                auto degree_map = right_core_degree_map->at(r);
                for(auto iter = degree_map->lower_bound(k); iter!=degree_map->end(); ++iter){
                    for(const auto&l:*iter->second){
                        if(!left_index_map->at(l)->count(k, j) || removed_l_set->count(l)){
                            continue;
                        }
                        if(!left_index_map->at(l)->count(k, j + 1) || previous_removed_l_set->count(l)){
                            if (!l_map->count(l)) {
                                auto l_degree = compute_left_vertex_core_degree(left_core_degree_map, right_index_map, l, k, j);
                                l_map->insert({l, l_degree});
                            }

                            --l_map->at(l);
                            if (l_map->at(l) < k) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
            removed_r_set->merge(*evicted_r_set);
        }
    }

    void branch_bipartite_core_maintenance::remove_left_vertex(const shared_ptr<abstract_bipartite_graph>& B,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_l_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_r_map,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                               const uint32_t l,
                                                               uint32_t i,
                                                               uint32_t j) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(l);

        while(!vertex_set->empty()){
            auto l1 = *vertex_set->begin();
            vertex_set->erase(l1);
            evicted_l_set->insert(l1);

            for(const auto &[r1, e1]:*B->get_left_vertex(l1)->get_edge_map()){
                if(candidate_r_map->count(r1)){
                    --candidate_r_map->at(r1);
                    if(candidate_r_map->at(r1) < j){
                        candidate_r_map->erase(r1);
                        evicted_r_set->insert(r1);

                        for(const auto &[l2,e2]:*B->get_right_vertex(r1)->get_edge_map()){
                            if(candidate_l_map->count(l2)){
                                --candidate_l_map->at(l2);
                                if(candidate_l_map->at(l2) < i){
                                    candidate_l_map->erase(l2);
                                    vertex_set->insert(l2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::remove_left_vertex(const shared_ptr<abstract_bipartite_graph>& B,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_l_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_r_map,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_l_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_r_map,
                                                               const uint32_t l,
                                                               uint32_t i,
                                                               uint32_t j) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(l);

        while(!vertex_set->empty()){
            auto l1 = *vertex_set->begin();
            vertex_set->erase(l1);
            evicted_l_set->insert(l1);

            for(const auto &[r1, e1]:*B->get_left_vertex(l1)->get_edge_map()){
                if(candidate_r_map->count(r1)){
                    --candidate_r_map->at(r1);
                    if(candidate_r_map->at(r1) < j){
                        candidate_r_map->erase(r1);
                        evicted_r_set->insert(r1);

                        for(const auto &[l2,e2]:*B->get_right_vertex(r1)->get_edge_map()){
                            if(candidate_l_map->count(l2)){
                                --candidate_l_map->at(l2);
                                if(candidate_l_map->at(l2) < i){
                                    candidate_l_map->erase(l2);
                                    vertex_set->insert(l2);
                                }
                            }
                            if(visited_l_map->count(l2)){
                                --visited_l_map->at(l2);
                            }
                        }
                    }
                }
                if(visited_r_map->count(r1)){
                    --visited_r_map->at(r1);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::remove_right_vertex(const shared_ptr<abstract_bipartite_graph>& B,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_l_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_r_map,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                                const uint32_t r,
                                                                uint32_t i,
                                                                uint32_t j)
    {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(r);

        while(!vertex_set->empty()){
            auto r1 = *vertex_set->begin();
            vertex_set->erase(r1);
            evicted_r_set->insert(r1);

            for(const auto &[l1,e1]:*B->get_right_vertex(r1)->get_edge_map()){
                if(candidate_l_map->count(l1)){
                    --candidate_l_map->at(l1);
                    if(candidate_l_map->at(l1) < i){
                        candidate_l_map->erase(l1);
                        evicted_l_set->insert(l1);

                        for(const auto &[r2, e2]:*B->get_left_vertex(l1)->get_edge_map()){
                            if(candidate_r_map->count(r2)){
                                --candidate_r_map->at(r2);
                                if(candidate_r_map->at(r2) < j){
                                    candidate_r_map->erase(r2);
                                    vertex_set->insert(r2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void branch_bipartite_core_maintenance::remove_right_vertex(const shared_ptr<abstract_bipartite_graph>& B,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_l_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_r_map,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_l_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> & visited_r_map,
                                                                const uint32_t r,
                                                                uint32_t i,
                                                                uint32_t j)
    {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(r);

        while(!vertex_set->empty()){
            auto r1 = *vertex_set->begin();
            vertex_set->erase(r1);
            evicted_r_set->insert(r1);

            for(const auto &[l1,e1]:*B->get_right_vertex(r1)->get_edge_map()){
                if(candidate_l_map->count(l1)){
                    --candidate_l_map->at(l1);
                    if(candidate_l_map->at(l1) < i){
                        candidate_l_map->erase(l1);
                        evicted_l_set->insert(l1);

                        for(const auto &[r2, e2]:*B->get_left_vertex(l1)->get_edge_map()){
                            if(candidate_r_map->count(r2)){
                                --candidate_r_map->at(r2);
                                if(candidate_r_map->at(r2) < j){
                                    candidate_r_map->erase(r2);
                                    vertex_set->insert(r2);
                                }
                            }

                            if(visited_r_map->count(r2)){
                                --visited_r_map->at(r2);
                            }
                        }
                    }
                }

                if(visited_l_map->count(l1)){
                    --visited_l_map->at(l1);
                }
            }
        }
    }


    void branch_bipartite_core_maintenance::remove_left_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_degree_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_degree_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                               uint32_t l,
                                                               uint32_t i,
                                                               uint32_t j) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(l);

        while(!vertex_set->empty()){
            auto l1 = *vertex_set->begin();
            vertex_set->erase(l1);
            evicted_l_set->insert(l1);

            auto l1_degree_map = left_degree_map->at(l1);
            for(auto iter1 = l1_degree_map->lower_bound(j); iter1 != l1_degree_map->end(); ++iter1){
                for(const auto &r1:*iter1->second){
                    if(inserted_r_map->count(r1)){
                        --inserted_r_map->at(r1);
                        if(inserted_r_map->at(r1) < j){
                            inserted_r_map->erase(r1);
                            evicted_r_set->insert(r1);

                            auto r1_degree_map = right_degree_map->at(r1);
                            for(auto iter2 = r1_degree_map->lower_bound(i);iter2!=r1_degree_map->end();++iter2){
                                for(const auto &l2:*iter2->second){
                                    if(inserted_l_map->count(l2)){
                                        --inserted_l_map->at(l2);
                                        if(inserted_l_map->at(l2) < i){
                                            inserted_l_map->erase(l2);
                                            vertex_set->insert(l2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::remove_vertex(const shared_ptr<abstract_bipartite_graph>& B,
                                                          uint32_t v,
                                                          const shared_ptr<unordered_set<uint32_t>>& evicted_l_set,
                                                          const shared_ptr<unordered_set<uint32_t>>& evicted_r_set,
                                                          const shared_ptr<unordered_set<uint32_t>>& inserted_l_set,
                                                          const shared_ptr<unordered_set<uint32_t>>& inserted_r_set,
                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& vertex_degree_map,
                                                          const shared_ptr<extend_list<int, uint32_t>>& order_list,
                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& k_rem_map,
                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>>& ext,
                                                          const shared_ptr<map<int, uint32_t>>& affected_vertex_map,
                                                          uint32_t i,
                                                          uint32_t j){
        auto pivot_node = order_list->find(v)->get_next();
        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                auto l_node = make_shared<extend_node<int, uint32_t>>(0, l);
                order_list->insert_before(l_node, pivot_node);
                k_rem_map->at(l) = vertex_degree_map->at(l);

                inserted_l_set->erase(l);
                vertex_degree_map->erase(l);

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (vertex_degree_map->count(r)) {
                        --vertex_degree_map->at(r);
                        if (vertex_degree_map->at(r) < j) {
                            evicted_r_set->insert(r);
                        }
                    }else if(ext->count(r)){
                        --ext->at(r);
                        if(ext->at(r) == 0){
                            ext->erase(r);
                            affected_vertex_map->erase(order_list->find_key(r).value());
                        }
                    }
                }
            }

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                auto r_node = make_shared<extend_node<int, uint32_t>>(0, r);
                order_list->insert_before(r_node, pivot_node);
                k_rem_map->at(r) = vertex_degree_map->at(r);

                inserted_r_set->erase(r);
                vertex_degree_map->erase(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (vertex_degree_map->count(l)) {
                        --vertex_degree_map->at(l);
                        if (vertex_degree_map->at(l) < i) {
                            evicted_l_set->insert(l);
                        }
                    }else if(ext->count(l)){
                        --ext->at(l);
                        if(ext->at(l) == 0){
                            ext->erase(l);
                            affected_vertex_map->erase(order_list->find_key(l).value());
                        }
                    }
                }
            }
        }
    }


    void branch_bipartite_core_maintenance::remove_left_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_degree_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_degree_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                               const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_l_map,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_r_map,
                                                               uint32_t l,
                                                               uint32_t i,
                                                               uint32_t j) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(l);

        while(!vertex_set->empty()){
            auto l1 = *vertex_set->begin();
            vertex_set->erase(l1);
            evicted_l_set->insert(l1);

            auto l1_degree_map = left_degree_map->at(l1);
            for(auto iter1 = l1_degree_map->lower_bound(j); iter1 != l1_degree_map->end(); ++iter1){
                for(const auto &r1:*iter1->second){
                    if(inserted_r_map->count(r1)){
                        --inserted_r_map->at(r1);
                        if(inserted_r_map->at(r1) < j){
                            inserted_r_map->erase(r1);
                            evicted_r_set->insert(r1);

                            auto r1_degree_map = right_degree_map->at(r1);
                            for(auto iter2 = r1_degree_map->lower_bound(i);iter2!=r1_degree_map->end();++iter2){
                                for(const auto &l2:*iter2->second){
                                    if(inserted_l_map->count(l2)){
                                        --inserted_l_map->at(l2);
                                        if(inserted_l_map->at(l2) < i){
                                            inserted_l_map->erase(l2);
                                            vertex_set->insert(l2);
                                        }
                                    }

                                    if(visited_l_map->count(l1)){
                                        --visited_l_map->at(l1);
                                    }
                                }
                            }
                        }
                    }
                    if(visited_r_map->count(r1)){
                        --visited_r_map->at(r1);
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::remove_right_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                                const uint32_t r,
                                                                uint32_t i,
                                                                uint32_t j)
    {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(r);

        while(!vertex_set->empty()){
            auto r1 = *vertex_set->begin();
            vertex_set->erase(r1);
            evicted_r_set->insert(r1);

            auto r1_degree_map = right_core_degree_map->at(r1);
            for(auto iter1 = r1_degree_map->lower_bound(i); iter1 != r1_degree_map->end(); ++iter1){
                for(const auto &l1:*iter1->second){
                    if(inserted_l_map->count(l1)){
                        --inserted_l_map->at(l1);
                        if(inserted_l_map->at(l1) < i){
                            inserted_l_map->erase(l1);
                            evicted_l_set->insert(l1);

                            auto l1_degree_map = left_core_degree_map->at(l1);
                            for(auto iter2 = l1_degree_map->lower_bound(j); iter2!=l1_degree_map->end();++iter2){
                                for(const auto&r2:*iter2->second){
                                    if(inserted_r_map->count(r2)){
                                        --inserted_r_map->at(r2);
                                        if(inserted_r_map->at(r2) < j){
                                            inserted_r_map->erase(r2);
                                            vertex_set->insert(r2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::remove_right_vertex(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_l_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& inserted_r_map,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                                const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_l_map,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>>& visited_r_map,
                                                                uint32_t r,
                                                                uint32_t i,
                                                                uint32_t j)
    {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(r);

        while(!vertex_set->empty()){
            auto r1 = *vertex_set->begin();
            vertex_set->erase(r1);
            evicted_r_set->insert(r1);

            auto r1_degree_map = right_core_degree_map->at(r1);
            for(auto iter1 = r1_degree_map->lower_bound(i); iter1 != r1_degree_map->end(); ++iter1){
                for(const auto &l1:*iter1->second){
                    if(inserted_l_map->count(l1)){
                        --inserted_l_map->at(l1);
                        if(inserted_l_map->at(l1) < i){
                            inserted_l_map->erase(l1);
                            evicted_l_set->insert(l1);

                            auto l1_degree_map = left_core_degree_map->at(l1);
                            for(auto iter2 = l1_degree_map->lower_bound(j); iter2!=l1_degree_map->end();++iter2){
                                for(const auto&r2:*iter2->second){
                                    if(inserted_r_map->count(r2)){
                                        --inserted_r_map->at(r2);
                                        if(inserted_r_map->at(r2) < j){
                                            inserted_r_map->erase(r2);
                                            vertex_set->insert(r2);
                                        }
                                    }

                                    if(visited_r_map->count(r2)){
                                        --visited_r_map->at(r2);
                                    }
                                }
                            }
                        }
                    }

                    if(visited_l_map->count(l1)){
                        --visited_l_map->at(l1);
                    }
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                           uint32_t max_i,
                                                                           uint32_t max_j){
        for (uint32_t i = max_i; i > 1; --i) {
            auto removed_l_set = make_shared<unordered_set<uint32_t>>();
            auto removed_r_set = make_shared<unordered_set<uint32_t>>();

            for (const auto &e: *edge_set) {
                if (left_index_map->at(e->get_left_vertex_id())->count(i, 1) &&
                    right_index_map->at(e->get_right_vertex_id())->count(1, i)) {
                    {
                        auto l1 = e->get_left_vertex_id();
                        if (B->get_left_vertex(l1) && B->get_left_vertex(l1)->get_degree() < i) {
                            removed_l_set->insert(l1);

                            for (const auto &[r1, e1]: *B->get_left_vertex(l1)->get_edge_map()) {
                                if (right_index_map->at(r1)->count(1, i)) {
                                    auto flag = false;
                                    for (const auto &[l2, e2]: *B->get_right_vertex(r1)->get_edge_map()) {
                                        if (B->get_left_vertex(l2)->get_degree() >= i) {
                                            flag = true;
                                        }
                                    }
                                    if (!flag) {
                                        removed_r_set->insert(r1);
                                    }
                                }
                            }
                        }
                    }

                    {
                        auto r1 = e->get_right_vertex_id();
                        if (B->get_right_vertex(r1)) {
                            auto flag = false;
                            for (const auto &[l1, e1]: *B->get_right_vertex(r1)->get_edge_map()) {
                                if (B->get_left_vertex(l1)->get_degree() >= i) {
                                    flag = true;
                                    break;
                                }
                            }
                            if (!flag) {
                                removed_r_set->insert(r1);
                            }
                        }
                    }
                }
            }

            for (const auto &l: *removed_l_set) {
                new_left_index_map->at(l)->remove(i, 1);
            }

            for (const auto &r: *removed_r_set) {
                new_right_index_map->at(r)->remove(1, i);
            }
        }
        for (uint32_t j = max_j; j > 1; --j) {
            auto removed_l_set = make_shared<unordered_set<uint32_t>>();
            auto removed_r_set = make_shared<unordered_set<uint32_t>>();

            for (const auto &e: *edge_set) {
                if (left_index_map->at(e->get_left_vertex_id())->count(1, j) &&
                    right_index_map->at(e->get_right_vertex_id())->count(j, 1)) {
                    {
                        auto r1 = e->get_right_vertex_id();
                        if (B->get_right_vertex(r1) && B->get_right_vertex(r1)->get_degree() < j) {
                            removed_r_set->insert(r1);

                            for (const auto &[l1, e1]: *B->get_right_vertex(r1)->get_edge_map()) {
                                if (left_index_map->at(l1)->count(1, j)) {
                                    auto flag = false;
                                    for (const auto &[r2, e2]: *B->get_left_vertex(l1)->get_edge_map()) {
                                        if (B->get_right_vertex(r2)->get_degree() >= j) {
                                            flag = true;
                                            break;
                                        }
                                    }
                                    if (!flag) {
                                        removed_l_set->insert(l1);
                                    }
                                }
                            }
                        }
                    }

                    {
                        auto l1 = e->get_left_vertex_id();
                        if (B->get_left_vertex(l1)) {
                            auto flag = false;
                            for (const auto &[r2, e2]: *B->get_left_vertex(l1)->get_edge_map()) {
                                if (B->get_right_vertex(r2)->get_degree() >= j) {
                                    flag = true;
                                }
                            }
                            if (!flag) {
                                removed_l_set->insert(l1);
                            }
                        }
                    }
                }
            }

            for (const auto &l: *removed_l_set) {
                new_left_index_map->at(l)->remove(1, j);
            }

            for (const auto &r: *removed_r_set) {
                new_right_index_map->at(r)->remove(j, 1);
            }
        }
    }

    void branch_bipartite_core_maintenance::update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                                           const shared_ptr<mutex>& global_left_mutex,
                                                                           const shared_ptr<mutex>& global_right_mutex,
                                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                           uint32_t max_i,
                                                                           uint32_t max_j,
                                                                           const shared_ptr<thread_pool>& pool){
        pool->submit_task([=] {
            for (uint32_t i = max_i; i > 1; --i) {
                auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto &e: *edge_set) {
                    if (left_index_map->at(e->get_left_vertex_id())->count(i, 1) &&
                        right_index_map->at(e->get_right_vertex_id())->count(1, i)) {
                        {
                            auto l1 = e->get_left_vertex_id();
                            if (B->get_left_vertex(l1) && B->get_left_vertex(l1)->get_degree() < i) {
                                removed_l_set->insert(l1);

                                for (const auto &[r1, e1]: *B->get_left_vertex(l1)->get_edge_map()) {
                                    if (right_index_map->at(r1)->count(1, i)) {
                                        auto flag = false;
                                        for (const auto &[l2, e2]: *B->get_right_vertex(r1)->get_edge_map()) {
                                            if (B->get_left_vertex(l2)->get_degree() >= i) {
                                                flag = true;
                                            }
                                        }
                                        if (!flag) {
                                            removed_r_set->insert(r1);
                                        }
                                    }
                                }
                            }
                        }

                        {
                            auto r1 = e->get_right_vertex_id();
                            if (B->get_right_vertex(r1)) {
                                auto flag = false;
                                for (const auto &[l1, e1]: *B->get_right_vertex(r1)->get_edge_map()) {
                                    if (B->get_left_vertex(l1)->get_degree() >= i) {
                                        flag = true;
                                        break;
                                    }
                                }
                                if (!flag) {
                                    removed_r_set->insert(r1);
                                }
                            }
                        }
                    }
                }

                global_left_mutex->lock();
                for (const auto &l: *removed_l_set) {
                    new_left_index_map->at(l)->remove(i, 1);
                }
                global_left_mutex->unlock();

                global_right_mutex->lock();
                for (const auto &r: *removed_r_set) {
                    new_right_index_map->at(r)->remove(1, i);
                }
                global_right_mutex->unlock();
            }
        });

        pool->submit_task([=] {
            for (uint32_t j = max_j; j > 1; --j) {
                auto removed_l_set = make_shared<unordered_set<uint32_t>>();
                auto removed_r_set = make_shared<unordered_set<uint32_t>>();

                for (const auto &e: *edge_set) {
                    if (left_index_map->at(e->get_left_vertex_id())->count(1, j) &&
                        right_index_map->at(e->get_right_vertex_id())->count(j, 1)) {
                        {
                            auto r1 = e->get_right_vertex_id();
                            if (B->get_right_vertex(r1) && B->get_right_vertex(r1)->get_degree() < j) {
                                removed_r_set->insert(r1);

                                for (const auto &[l1, e1]: *B->get_right_vertex(r1)->get_edge_map()) {
                                    if (left_index_map->at(l1)->count(1, j)) {
                                        auto flag = false;
                                        for (const auto &[r2, e2]: *B->get_left_vertex(l1)->get_edge_map()) {
                                            if (B->get_right_vertex(r2)->get_degree() >= j) {
                                                flag = true;
                                                break;
                                            }
                                        }
                                        if (!flag) {
                                            removed_l_set->insert(l1);
                                        }
                                    }
                                }
                            }
                        }

                        {
                            auto l1 = e->get_left_vertex_id();
                            if (B->get_left_vertex(l1)) {
                                auto flag = false;
                                for (const auto &[r2, e2]: *B->get_left_vertex(l1)->get_edge_map()) {
                                    if (B->get_right_vertex(r2)->get_degree() >= j) {
                                        flag = true;
                                    }
                                }
                                if (!flag) {
                                    removed_l_set->insert(l1);
                                }
                            }
                        }
                    }
                }

                global_left_mutex->lock();
                for (const auto &l: *removed_l_set) {
                    new_left_index_map->at(l)->remove(1, j);
                }
                global_left_mutex->unlock();

                global_right_mutex->lock();
                for (const auto &r: *removed_r_set) {
                    new_right_index_map->at(r)->remove(j, 1);
                }
                global_right_mutex->unlock();
            }
        });
        pool->barrier();
    }

    void branch_bipartite_core_maintenance::update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map){
        auto affected_l_set = make_shared<unordered_set<uint64_t>>();
        auto affected_r_set = make_shared<unordered_set<uint64_t>>();

        auto inserted_l_set = make_shared<unordered_set<uint64_t>>();
        auto inserted_r_set = make_shared<unordered_set<uint64_t>>();

        B->insert_edge_collection(edge_set);

        for(const auto&e:*edge_set){
            auto l = e->get_left_vertex_id();
            left_index_map->at(l)->remove(B->get_left_vertex(l)->get_degree());
            if(!affected_l_set->count(l)){
                affected_l_set->insert(l);
                inserted_l_set->insert(l);
                for(const auto &[r1,e1]:*B->get_left_vertex(l)->get_edge_map()){
                    inserted_r_set->insert(r1);
                }
            }

            auto r = e->get_right_vertex_id();
            right_index_map->at(r)->remove(B->get_right_vertex(r)->get_degree());
            if(!affected_r_set->count(r)){
                affected_r_set->insert(r);
                inserted_r_set->insert(r);
                for(const auto &[l1, e1]:*B->get_right_vertex(r)->get_edge_map()){
                    inserted_l_set->insert(l1);
                }
            }
            B->remove_edge(e);
        }

        for(const auto &l:*inserted_l_set){
            if(!B->get_left_vertex(l)){
                continue;
            }
            uint64_t max_degree = 0;
            for(const auto &[r, e]:*B->get_left_vertex(l)->get_edge_map()){
                if(B->get_right_vertex(r)->get_degree() > max_degree){
                    max_degree = B->get_right_vertex(r)->get_degree();
                }
            }
            left_index_map->at(l)->remove(1, max_degree);
        }

        for(const auto &r:*inserted_r_set){
            if(!B->get_right_vertex(r)){
                continue;
            }
            uint64_t max_degree = 0;
            for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                if(B->get_left_vertex(l)->get_degree() > max_degree){
                    max_degree = B->get_left_vertex(l)->get_degree();
                }
            }
            right_index_map->at(r)->remove(1, max_degree);
        }
    }


    void branch_bipartite_core_maintenance::update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                           const shared_ptr<thread_pool>& pool){
        auto affected_l_set = make_shared<unordered_set<uint64_t>>();
        auto affected_r_set = make_shared<unordered_set<uint64_t>>();

        auto inserted_l_set = make_shared<unordered_set<uint64_t>>();
        auto inserted_r_set = make_shared<unordered_set<uint64_t>>();

        B->insert_edge_collection(edge_set);

        for(const auto&e:*edge_set){
            auto l = e->get_left_vertex_id();
            left_index_map->at(l)->remove(B->get_left_vertex(l)->get_degree());
            if(!affected_l_set->count(l)){
                affected_l_set->insert(l);
                inserted_l_set->insert(l);
                for(const auto &[r1,e1]:*B->get_left_vertex(l)->get_edge_map()){
                    inserted_r_set->insert(r1);
                }
            }

            auto r = e->get_right_vertex_id();
            right_index_map->at(r)->remove(B->get_right_vertex(r)->get_degree());
            if(!affected_r_set->count(r)){
                affected_r_set->insert(r);
                inserted_r_set->insert(r);
                for(const auto &[l1, e1]:*B->get_right_vertex(r)->get_edge_map()){
                    inserted_l_set->insert(l1);
                }
            }
            B->remove_edge(e);
        }

        for(const auto &l:*inserted_l_set){
            if(!B->get_left_vertex(l)){
                continue;
            }
            pool->submit_task([=]{
                uint64_t max_degree = 0;
                for(const auto &[r, e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(B->get_right_vertex(r)->get_degree() > max_degree){
                        max_degree = B->get_right_vertex(r)->get_degree();
                    }
                }
                left_index_map->at(l)->remove(1, max_degree);
            });
        }

        for(const auto &r:*inserted_r_set){
            if(!B->get_right_vertex(r)){
                continue;
            }
            pool->submit_task([=]{
                uint64_t max_degree = 0;
                for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(B->get_left_vertex(l)->get_degree() > max_degree){
                        max_degree = B->get_left_vertex(l)->get_degree();
                    }
                }
                right_index_map->at(r)->remove(1, max_degree);
            });
        }
        pool->barrier();
    }


    void branch_bipartite_core_maintenance::update_index_for_insertion(const shared_ptr<abstract_bipartite_graph> &B,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                       const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                                       const shared_ptr<unordered_set<uint32_t>> &inserted_r_set){
        /**
         * @brief remove update
         */
        B->remove_edge_collection(edge_set);
        for(const auto &l:*inserted_l_set) {
            auto l_vertex = B->get_left_vertex(l);
            if(l_vertex){
                auto degree = l_vertex->get_degree();
                for (const auto &[r, e]: *l_vertex->get_edge_map()) {

                    auto r_degree_map = right_core_degree_map->at(r);
                    r_degree_map->at(degree)->erase(l);
                    if (r_degree_map->at(degree)->empty()) {
                        r_degree_map->erase(degree);
                    }
                }
            }else{
                left_core_degree_map->insert({l, make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
            }
        }

        for(const auto &r:*inserted_r_set){
            auto r_vertex = B->get_right_vertex(r);
            if(r_vertex){
                auto degree = r_vertex->get_degree();
                for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                    auto l_degree_map = left_core_degree_map->at(l);
                    l_degree_map->at(degree)->erase(r);
                    if (l_degree_map->at(degree)->empty()) {
                        l_degree_map->erase(degree);
                    }
                }
            }else{
                right_core_degree_map->insert({r, make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
            }
        }

        B->insert_edge_collection(edge_set);
        /**
         * @brief insert update
         */
        for(const auto &l:*inserted_l_set) {
            auto l_vertex = B->get_left_vertex(l);
            auto degree = l_vertex->get_degree();
            for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                auto r_degree_map = right_core_degree_map->at(r);
                if (!r_degree_map->count(degree)) {
                    r_degree_map->insert({degree, make_shared<unordered_set<uint32_t>>()});
                }
                r_degree_map->at(degree)->insert(l);
            }
        }

        for(const auto &r:*inserted_r_set){
            auto r_vertex = B->get_right_vertex(r);
            auto degree = r_vertex->get_degree();
            for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                auto l_degree_map = left_core_degree_map->at(l);
                if (!l_degree_map->count(degree)) {
                    l_degree_map->insert({degree, make_shared<unordered_set<uint32_t>>()});
                }
                l_degree_map->at(degree)->insert(r);
            }
        }
    }


    void branch_bipartite_core_maintenance::update_index_for_removal(const shared_ptr<abstract_bipartite_graph> &B,
                                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                     const shared_ptr<unordered_set<uint32_t>> &removed_l_set,
                                                                     const shared_ptr<unordered_set<uint32_t>> &removed_r_set){

        /**
        * @brief find removed vertices and reinsert vertices
        */
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        B->insert_edge_collection(edge_set);
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            l_set->insert(l);
            r_set->insert(r);
        }

        /**
         * @brief remove update
         */
        for(const auto &l:*l_set) {
            auto l_vertex = B->get_left_vertex(l);
            auto degree = l_vertex->get_degree();
            for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                if(!edge_set->count(e)){
                    auto r_degree_map = right_core_degree_map->at(r);
                    r_degree_map->at(degree)->erase(l);
                    if (r_degree_map->at(degree)->empty()) {
                        r_degree_map->erase(degree);
                    }
                }
            }
        }

        for(const auto &r:*r_set){
            auto r_vertex = B->get_right_vertex(r);
            auto degree = r_vertex->get_degree();
            for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                if(!edge_set->count(e)){
                    auto l_degree_map = left_core_degree_map->at(l);
                    l_degree_map->at(degree)->erase(r);
                    if (l_degree_map->at(degree)->empty()) {
                        l_degree_map->erase(degree);
                    }
                }
            }
        }

        /**
         * @brief insert update
         */

        B->remove_edge_collection(edge_set);


        for(const auto &l:*l_set){
            if (removed_l_set->count(l)) {
                continue;
            }
            auto l_vertex = B->get_left_vertex(l);
            auto l_degree = l_vertex->get_degree();
            for (const auto &[r, e]: *l_vertex->get_edge_map()) {
                if (!removed_r_set->count(r)) {

                    auto r_degree_map = right_core_degree_map->at(r);
                    if (!r_degree_map->count(l_degree)) {
                        r_degree_map->insert({l_degree, make_shared<unordered_set<uint32_t>>()});
                    }
                    r_degree_map->at(l_degree)->insert(l);
                }
            }
        }

        for(const auto &r:*r_set){
            if (removed_r_set->count(r)) {
                continue;
            }
            auto r_vertex = B->get_right_vertex(r);
            auto r_degree = r_vertex->get_degree();
            for (const auto &[l, e]: *r_vertex->get_edge_map()) {
                if (!removed_l_set->count(l)) {

                    auto l_degree_map = left_core_degree_map->at(l);
                    if (!l_degree_map->count(r_degree)) {
                        l_degree_map->insert({r_degree, make_shared<unordered_set<uint32_t>>()});
                    }
                    l_degree_map->at(r_degree)->insert(r);
                }
            }
        }
    }



    void branch_bipartite_core_maintenance::update_left_core_degree_map(
            const shared_ptr<scnu::abstract_bipartite_graph> &B,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_left_store_index>>> &left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_right_store_index>>> &right_index_map,
            const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
            const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &core_degree_map,
            uint32_t i, uint32_t k) {
        for (const auto &l: *inserted_l_set) {
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (right_index_map->at(r)->count(k, i) && !right_index_map->at(r)->count(k, i + 1)) {
                    ++core_degree_map->at(r);
                }
            }
        }

        for (const auto &r: *inserted_r_set) {
            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if (left_index_map->at(l)->count(i, k) && !left_index_map->at(l)->count(i + 1, k)){
                    ++core_degree_map->at(l);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::update_middle_core_degree_map(
            const shared_ptr<scnu::abstract_bipartite_graph> &B,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_left_store_index>>> &left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_right_store_index>>> &right_index_map,
            const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
            const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &core_degree_map,
            uint32_t k) {
        for (const auto &l: *inserted_l_set) {
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (right_index_map->at(r)->count(k, k) && !right_index_map->at(r)->count(k + 1, k + 1)) {
                    ++core_degree_map->at(r);
                }
            }
        }

        for(const auto &r:*inserted_r_set){
            for(const auto &[l, e]:*B->get_right_vertex(r)->get_edge_map()) {
                if (left_index_map->at(l)->count(k, k) && !left_index_map->at(l)->count(k + 1, k + 1)) {
                    ++core_degree_map->at(l);
                }
            }
        }
    }

    void
    branch_bipartite_core_maintenance::update_right_core_degree_map(const shared_ptr<scnu::abstract_bipartite_graph> &B,
                                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_left_store_index>>> &left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<scnu::bipartite_core_right_store_index>>> &right_index_map,
                                                                    const shared_ptr<unordered_set<uint32_t>> &inserted_l_set,
                                                                    const shared_ptr<unordered_set<uint32_t>> &inserted_r_set,
                                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &core_degree_map,
                                                                    uint32_t k, uint32_t j) {
        for (const auto &l: *inserted_l_set) {
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (right_index_map->at(r)->count(j, k) && !right_index_map->at(r)->count(j + 1, k)) {
                    ++core_degree_map->at(r);
                }
            }
        }

        for (const auto &r: *inserted_r_set) {
            for(const auto &[l, e]:*B->get_right_vertex(r)->get_edge_map()){
                if(left_index_map->at(l)->count(k, j) && !left_index_map->at(l)->count(k, j + 1)){
                    ++core_degree_map->at(l);
                }
            }
        }
    }
}

