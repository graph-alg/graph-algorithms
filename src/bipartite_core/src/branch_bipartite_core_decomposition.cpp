
#include "bipartite_core/branch_bipartite_core_decomposition.h"

namespace scnu {

    uint32_t branch_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map) {
        auto left_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            auto degree = l_vertex->get_degree();
            left_map->insert({l, degree});

            left_index_map->insert({l, make_shared<left_vertex_index>()});
            for(uint32_t i = 1; i<= degree; ++i){
                left_index_map->at(l)->insert(i, 1);
            }
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            auto degree = r_vertex->get_degree();
            right_map->insert({r, degree});

            right_index_map->insert({r, make_shared<right_vertex_index>()});
            for(uint32_t j = 1; j <= degree; ++j){
                right_index_map->at(r)->insert(j, 1);
            }
        }

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            auto degree = l_vertex->get_degree();
            for(const auto &[r, e]:*l_vertex->get_edge_map()){
                right_index_map->at(r)->insert(1, degree);
            }
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            auto degree = r_vertex->get_degree();
            for(const auto &[l, e]:*r_vertex->get_edge_map()){
                left_index_map->at(l)->insert(1, degree);
            }
        }


        uint32_t max_k = 1;

        for (uint32_t k = 2; true; ++k) {
            {
                find_core(B, left_map, right_map, k);
                if(!left_map->empty() && !right_map->empty()) {
                    max_k = k;
                }else{
                    break;
                }
            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                find_left_core(B, sub_left_degree_map, sub_right_degree_map,
                                   left_index_map,
                                   right_index_map, k);
            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                find_right_core(B, sub_left_degree_map, sub_right_degree_map, left_index_map,
                                    right_index_map,
                                    k);
            }
        }

        return max_k;
    }

    uint32_t branch_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &B,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                                            const shared_ptr<thread_pool>& pool) {
        auto left_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto left_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto right_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            auto degree = l_vertex->get_degree();
            left_map->insert({l, degree});

            left_mutex_map->insert({l, shared_ptr<mutex>()});
            left_index_map->insert({l, shared_ptr<left_vertex_index>()});
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            auto degree = r_vertex->get_degree();
            right_map->insert({r, degree});

            right_mutex_map->insert({r, shared_ptr<mutex>()});
            right_index_map->insert({r, shared_ptr<right_vertex_index>()});
        }

        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=]{
                left_mutex_map->at(l) =  make_shared<mutex>();
                left_index_map->at(l) =  make_shared<left_vertex_index>();

                auto degree = l_vertex->get_degree();
                for(uint32_t i = 1; i<= degree; ++i){
                    left_index_map->at(l)->insert(i, 1);
                }
            });
        }

        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=]{
                right_mutex_map->at(r) =  make_shared<mutex>();
                right_index_map->at(r) = make_shared<right_vertex_index>();

                auto degree =  r_vertex->get_degree();
                for(uint32_t j = 1; j <= degree; ++j){
                    right_index_map->at(r)->insert(j, 1);
                }
            });
        }
        pool->barrier();
        for (const auto&[l, l_vertex]: *B->get_left_vertex_map()) {
            pool->submit_task([=]{
                auto degree = l_vertex->get_degree();
                for(const auto &[r, e]:*l_vertex->get_edge_map()){
                    right_mutex_map->at(r)->lock();
                    right_index_map->at(r)->insert(1, degree);
                    right_mutex_map->at(r)->unlock();
                }
            });
        }
        for (const auto&[r, r_vertex]: *B->get_right_vertex_map()) {
            pool->submit_task([=]{
                auto degree =  r_vertex->get_degree();
                for(const auto &[l, e]:*r_vertex->get_edge_map()){
                    left_mutex_map->at(l)->lock();
                    left_index_map->at(l)->insert(1, degree);
                    left_mutex_map->at(l)->unlock();
                }
            });
        }

        uint32_t max_k = 1;

        for (uint32_t k = 2; true; ++k) {
            {
                find_core(B, left_map, right_map, k);
                if(!left_map->empty() && !right_map->empty()) {
                    max_k = k;
                }else{
                    break;
                }
            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                pool->submit_task([=] {
                    find_left_core(B, left_mutex_map, right_mutex_map, sub_left_degree_map, sub_right_degree_map,
                                   left_index_map,
                                   right_index_map, k);
                });

            }

            {
                auto sub_left_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(left_map);
                auto sub_right_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(right_map);
                pool->submit_task([=] {
                    find_right_core(B, left_mutex_map, right_mutex_map, sub_left_degree_map, sub_right_degree_map, left_index_map,
                                    right_index_map,
                                    k);
                });
            }

        }
        pool->barrier();

        return max_k;
    }

    void branch_bipartite_core_decomposition::find_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_degree_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_degree_map,
                                                        uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();


        for (const auto&[l, l_degree]: *left_degree_map) {
            if (l_degree < k) {
                evicted_l_set->emplace(l);
            }
        }

        for (const auto &[r, r_degree]: *right_degree_map) {
            if (r_degree < k) {
                evicted_r_set->emplace(r);
            }
        }


        while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_degree_map->erase(l);

                for (const auto&[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_degree_map->count(r) && right_degree_map->at(r) >= k) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) < k) {
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }


            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_degree_map->erase(r);

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_degree_map->count(l) && left_degree_map->at(l) >= k) {
                        --left_degree_map->at(l);
                        if (left_degree_map->at(l) < k) {
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }
    }

    uint32_t branch_bipartite_core_decomposition::find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>>& right_index_map,
                                                                 uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_i = k;
        while (!left_map->empty()) {
            uint32_t  i = left_map->begin()->second;

            for (const auto&[l, l_degree]: *left_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == i)
                {
                    evicted_l_set->insert(l);
                }
            }

            max_i = i;

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_map->erase(l);

                for (uint32_t index = k; index <= i; ++index) {
                    left_index_map->at(l)->insert(index, k);
                }

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_map->count(r) && right_map->at(r) >= k) {
                        --right_map->at(r);
                        if (right_map->at(r) < k) {
                            right_map->erase(r);
                            right_index_map->at(r)->insert(k, i);

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (left_map->count(l2) && left_map->at(l2) > i) {
                                    --left_map->at(l2);
                                    if (left_map->at(l2) == i) {
                                        evicted_l_set->emplace(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_i;
    }


    uint32_t branch_bipartite_core_decomposition::find_left_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& left_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& right_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>>& right_index_map,
                                                                 uint32_t k) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_i = k;
        while (!left_map->empty()) {
            uint32_t  i = left_map->begin()->second;

            for (const auto&[l, l_degree]: *left_map) {
                if (l_degree < i) {
                    i = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == i)
                {
                    evicted_l_set->insert(l);
                }
            }

            max_i = i;

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_map->erase(l);

                left_mutex_map->at(l)->lock();
                for (uint32_t index = k; index <= i; ++index) {
                    left_index_map->at(l)->insert(index, k);
                }
                left_mutex_map->at(l)->unlock();

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_map->count(r) && right_map->at(r) >= k) {
                        --right_map->at(r);
                        if (right_map->at(r) < k) {
                            right_map->erase(r);

                            right_mutex_map->at(r)->lock();
                            right_index_map->at(r)->insert(k, i);
                            right_mutex_map->at(r)->unlock();

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (left_map->count(l2) && left_map->at(l2) > i) {
                                    --left_map->at(l2);
                                    if (left_map->at(l2) == i) {
                                        evicted_l_set->emplace(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_i;
    }

    uint32_t branch_bipartite_core_decomposition::find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                  const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>>& left_index_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>>& right_index_map,
                                                                  uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_j = k;
        while (!right_map->empty()) {
            uint32_t  j = right_map->begin()->second;

            for (const auto&[r, r_degree]: *right_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if(r_degree == j)
                {
                    evicted_r_set->insert(r);
                }
            }
            max_j = j;

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_map->erase(r);

                for (uint32_t index = k; index <= j; ++index) {
                    right_index_map->at(r)->insert(index, k);
                }

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_map->count(l) && left_map->at(l) >= k) {
                        --left_map->at(l);
                        if (left_map->at(l) < k) {
                            left_map->erase(l);
                            left_index_map->at(l)->insert(k, j);

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (right_map->count(r2) && right_map->at(r2) > j) {
                                    --right_map->at(r2);
                                    if (right_map->at(r2) == j) {
                                        evicted_r_set->insert(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_j;
    }


    uint32_t branch_bipartite_core_decomposition::find_right_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& left_mutex_map,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& right_mutex_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_map,
                                                                 const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>>& left_index_map,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>>& right_index_map,
                                                                 uint32_t k) {
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_j = k;
        while (!right_map->empty()) {
            uint32_t  j = right_map->begin()->second;

            for (const auto&[r, r_degree]: *right_map) {
                if (r_degree < j) {
                    j = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                } else if(r_degree == j)
                {
                    evicted_r_set->insert(r);
                }
            }
            max_j = j;

            while (!evicted_r_set->empty()) {
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);

                right_map->erase(r);

                right_mutex_map->at(r)->lock();
                for (uint32_t index = k; index <= j; ++index) {
                    right_index_map->at(r)->insert(index, k);
                }
                right_mutex_map->at(r)->unlock();

                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (left_map->count(l) && left_map->at(l) >= k) {
                        --left_map->at(l);
                        if (left_map->at(l) < k) {
                            left_mutex_map->at(l)->lock();
                            left_index_map->at(l)->insert(k, j);
                            left_mutex_map->at(l)->unlock();

                            for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                                if (right_map->count(r2) && right_map->at(r2) > j) {
                                    --right_map->at(r2);
                                    if (right_map->at(r2) == j) {
                                        evicted_r_set->emplace(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return max_j;
    }
}