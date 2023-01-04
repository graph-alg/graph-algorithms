
#include "bipartite_core/branch_bipartite_core_maintenance.h"

namespace scnu{

    void branch_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, make_shared<left_vertex_index>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()) {
            new_right_index_map->insert({r, make_shared<right_vertex_index>()});
        }
    }

    void branch_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
                                                 const shared_ptr<thread_pool>& pool) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, shared_ptr<left_vertex_index>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()){
            new_right_index_map->insert({r, shared_ptr<right_vertex_index>()});
        }

        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            pool->submit_task([=]{
                new_left_index_map->at(l) = make_shared<left_vertex_index>();
            });
        }
        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            pool->submit_task([=]{
                new_right_index_map->at(r) = make_shared<right_vertex_index>();
            });
        }
        pool->barrier();
    }


    void branch_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
                                                 const shared_ptr<thread_pool>& pool) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, shared_ptr<left_vertex_index>()});
            left_core_degree_map->insert({l, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()){
            new_right_index_map->insert({r, shared_ptr<right_vertex_index>()});
            right_core_degree_map->insert({r, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }

        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            pool->submit_task([=]{
                new_left_index_map->at(l) = make_shared<left_vertex_index>();
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
                new_right_index_map->at(r) = make_shared<right_vertex_index>();
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
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& new_right_index_map){

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
                left_index_map->insert({l, make_shared<left_vertex_index>()});
                new_left_index_map->insert({l, make_shared<left_vertex_index>()});

                previous_middle_inserted_l_map->insert({l, 0});
            }
            left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);

            affected_r_set->insert(r);
            if (!right_index_map->count(r)) {
                right_index_map->insert({r, make_shared<right_vertex_index>()});
                new_right_index_map->insert({r, make_shared<right_vertex_index>()});
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
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& new_right_index_map,
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
                    left_index_map->insert({l, make_shared<left_vertex_index>()});
                    new_left_index_map->insert({l, make_shared<left_vertex_index>()});
                }
                left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);


                if (!right_index_map->count(r)) {
                    right_index_map->insert({r, make_shared<right_vertex_index>()});
                    new_right_index_map->insert({r, make_shared<right_vertex_index>()});
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

        for(const auto &[r, r_index]:*new_right_index_map){
            for(const auto &[j, i]:*r_index->get_index_map()){
                right_index_map->at(r)->insert(j, i);
            }
            r_index->clear();
        }
    }


    void branch_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& new_right_index_map,
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
                    left_index_map->insert({l, make_shared<left_vertex_index>()});
                    new_left_index_map->insert({l, make_shared<left_vertex_index>()});
                }
                left_index_map->at(l)->insert(B->get_left_vertex(l)->get_degree(), 1);


                if (!right_index_map->count(r)) {
                    right_index_map->insert({r, make_shared<right_vertex_index>()});
                    new_right_index_map->insert({r, make_shared<right_vertex_index>()});
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

        for(const auto &[r, r_index]:*new_right_index_map){
            for(const auto &[j, i]:*r_index->get_index_map()){
                 right_index_map->at(r)->insert(j, i);
            }
            r_index->clear();
        }
    }

    void branch_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map) {
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
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                    const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
                                                   const shared_ptr<thread_pool>& pool) {
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
                                                            const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
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


    uint32_t branch_bipartite_core_maintenance::compute_left_vertex_core_degree(const shared_ptr<abstract_bipartite_graph>&B,
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_vertex_index_map,
                                                                                uint32_t l,
                                                                                uint32_t i,
                                                                                uint32_t j)
    {
        uint32_t core_degree = 0;
        for(const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()){
            if(right_vertex_index_map->at(r)->count(j, i)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph>&B,
                                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_vertex_index_map,
                                                                                 uint32_t r,
                                                                                 uint32_t i,
                                                                                 uint32_t j){
        uint32_t core_degree = 0;
        for(const auto&[l,e]:*B->get_right_vertex(r)->get_edge_map()){
            if(left_vertex_index_map->at(l)->count(i, j)){
                ++core_degree;
            }
        }
        return core_degree;
    }

    uint32_t branch_bipartite_core_maintenance::compute_left_vertex_core_degree(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
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
                                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
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
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_vertex_index_map,
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
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_vertex_index_map,
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
                                                                                const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_vertex_index_map,
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
                                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>>& left_vertex_index_map,
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
                                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>>& left_vertex_index_map,
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
                                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>>& left_vertex_index_map,
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
                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
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
                } else
                {
                    remove_right_vertex(B, inserted_l_map, inserted_r_map,
                                        evicted_l_set, evicted_r_set, visited_l_map, visited_r_map, r, i, k);
                }
            }
        }
    }


    void branch_bipartite_core_maintenance::left_candidate_graph(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
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
                                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
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
                if(degree >= k){
                    inserted_r_map->insert({r,degree});
                } else
                {
                    remove_right_vertex(B, inserted_l_map, inserted_r_map,evicted_l_set, evicted_r_set,
                                        visited_l_set, visited_r_set, r, k, k);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::middle_candidate_graph(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_l_map,
                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &previous_inserted_r_map,
                                                                   uint32_t k,
                                                                   const shared_ptr<unordered_map<uint32_t,uint32_t>> &inserted_l_map,
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
                                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map,
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
                if(degree >= j){
                    inserted_r_map->insert({r, degree});
                } else
                {
                    remove_right_vertex(B, inserted_l_map, inserted_r_map, evicted_l_set, evicted_r_set,
                                        visited_l_map, visited_r_map, r, k, j);
                }
            }
        }
    }

    void branch_bipartite_core_maintenance::right_candidate_graph(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_vertex_index_map,
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
                                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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

    void branch_bipartite_core_maintenance::left_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                            const shared_ptr<mutex>& global_left_mutex,
                                                                            const shared_ptr<mutex>& global_right_mutex,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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


    void branch_bipartite_core_maintenance::middle_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                              const shared_ptr<mutex>& global_left_mutex,
                                                                              const shared_ptr<mutex>& global_right_mutex,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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

    void branch_bipartite_core_maintenance::right_partial_core_decomposition(const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                             const shared_ptr<mutex>& global_left_mutex,
                                                                             const shared_ptr<mutex>& global_right_mutex,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_l_map,
                                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &inserted_r_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
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

    void branch_bipartite_core_maintenance::left_removal_partial_core(const shared_ptr<abstract_bipartite_graph>& B,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
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
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
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

    void branch_bipartite_core_maintenance::middle_removal_partial_core(const shared_ptr<abstract_bipartite_graph>& B,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
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
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
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

    void branch_bipartite_core_maintenance::right_removal_partial_core(const shared_ptr<abstract_bipartite_graph>& B,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &left_core_degree_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &right_core_degree_map,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
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
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &new_left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &new_right_index_map,
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
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map){
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
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
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
}

