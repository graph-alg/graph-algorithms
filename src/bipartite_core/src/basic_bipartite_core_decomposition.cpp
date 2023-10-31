
#include "bipartite_core/basic_bipartite_core_decomposition.h"

namespace scnu {

    shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>
    basic_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map) {

        auto core_left_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        auto core_right_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();

        core_left_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
        core_right_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
        for (const auto&[l, l_vertex]: *G->get_left_vertex_map()) {
            core_left_vertex_degree_map->at({1, 1})->insert({l, l_vertex->get_degree()});
            left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
            left_index_map->at(l)->insert(1, 1);
        }
        for (const auto&[r, r_vertex]: *G->get_right_vertex_map()) {
            core_right_vertex_degree_map->at({1, 1})->insert({r, r_vertex->get_degree()});
            right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
            right_index_map->at(r)->insert(1, 1);
        }

        auto total_pair_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
        total_pair_set->insert({1, 1});

        auto current_pair_queue = make_shared<queue<pair<uint32_t, uint32_t>>>();
        current_pair_queue->push({1,1});

        while (!current_pair_queue->empty()) {
            auto p  = current_pair_queue->front();
            current_pair_queue->pop();

            auto left_vertex_degree_map = core_left_vertex_degree_map->at(p);
            auto right_vertex_degree_map = core_right_vertex_degree_map->at(p);

            auto[i, j] = p;

            if (!total_pair_set->count({i + 1, j}) && !total_pair_set->count({i, j + 1})) {
                {
                    auto sub_left_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                            left_vertex_degree_map);
                    auto sub_right_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                            right_vertex_degree_map);

                    basic_bipartite_core_decomposition::find_left_core(G, sub_left_vertex_degree_map,
                                                                       sub_right_vertex_degree_map, i + 1, j);
                    if (!sub_left_vertex_degree_map->empty()) {
                        for (const auto&[l, l_degree]: *sub_left_vertex_degree_map) {
                            left_index_map->at(l)->insert(i + 1, j);
                        }
                        for (const auto&[r, r_degree]: *sub_right_vertex_degree_map) {
                            right_index_map->at(r)->insert(j, i + 1);
                        }

                        total_pair_set->insert({i + 1, j});
                        current_pair_queue->push({i + 1, j});

                        core_left_vertex_degree_map->insert({make_pair(i + 1, j), sub_left_vertex_degree_map});
                        core_right_vertex_degree_map->insert({make_pair(i + 1, j), sub_right_vertex_degree_map});
                    }
                }
                {
                    find_right_core(G, left_vertex_degree_map, right_vertex_degree_map, i, j + 1);
                    if (!left_vertex_degree_map->empty()) {
                        for (const auto&[l, l_degree]: *left_vertex_degree_map) {
                            left_index_map->at(l)->insert(i, j + 1);
                        }
                        for (const auto&[r, r_degree]: *right_vertex_degree_map) {
                            right_index_map->at(r)->insert(j + 1, i);
                        }


                        total_pair_set->insert({i, j + 1});
                        current_pair_queue->push({i, j + 1});

                        core_left_vertex_degree_map->insert({make_pair(i, j + 1), left_vertex_degree_map});
                        core_right_vertex_degree_map->insert({make_pair(i, j + 1), right_vertex_degree_map});
                    }
                }
            } else if (!total_pair_set->count({i + 1, j})) {
                find_left_core(G, left_vertex_degree_map, right_vertex_degree_map, i + 1, j);
                if (!left_vertex_degree_map->empty()) {
                    for (const auto&[l, l_degree]: *left_vertex_degree_map) {
                        left_index_map->at(l)->insert(i + 1, j);
                    }
                    for (const auto&[r, r_degree]: *right_vertex_degree_map) {
                        right_index_map->at(r)->insert(j, i + 1);
                    }

                    total_pair_set->insert({i + 1, j});
                    current_pair_queue->push({i + 1, j});

                    core_left_vertex_degree_map->insert({{i + 1, j}, left_vertex_degree_map});
                    core_right_vertex_degree_map->insert({{i + 1, j}, right_vertex_degree_map});
                }
            } else if (!total_pair_set->count({i, j + 1})) {
                find_right_core(G, left_vertex_degree_map, right_vertex_degree_map, i, j + 1);
                if (!left_vertex_degree_map->empty()) {
                    for (const auto&[l, l_degree]: *left_vertex_degree_map) {
                        left_index_map->at(l)->insert(i, j + 1);
                    }
                    for (const auto&[r, r_degree]: *right_vertex_degree_map) {
                        right_index_map->at(r)->insert(j + 1, i);
                    }

                    total_pair_set->insert({i, j + 1});
                    current_pair_queue->push({i, j + 1});

                    core_left_vertex_degree_map->insert({make_pair(i, j + 1), left_vertex_degree_map});
                    core_right_vertex_degree_map->insert({make_pair(i, j + 1), right_vertex_degree_map});
                }
            }
        }

        return total_pair_set;
    }


    shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>
    basic_bipartite_core_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                  uint32_t thread_number) {

        auto core_mutex = make_shared<mutex>();
        auto core_left_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        auto core_right_vertex_degree_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();

        auto left_vertex_mutex_map = make_shared<unordered_map<uint32_t,shared_ptr<mutex>>>();
        auto right_vertex_mutex_map = make_shared<unordered_map<uint32_t,shared_ptr<mutex>>>();

        core_left_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
        core_right_vertex_degree_map->insert({{1, 1}, make_shared<unordered_map<uint32_t, uint32_t>>()});
        for (const auto&[l, l_vertex]: *G->get_left_vertex_map()) {
            core_left_vertex_degree_map->at({1, 1})->insert({l, l_vertex->get_degree()});
            left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
            left_index_map->at(l)->insert(1, 1);

            left_vertex_mutex_map->insert({l, make_shared<mutex>()});
        }
        for (const auto&[r, r_vertex]: *G->get_right_vertex_map()) {
            core_right_vertex_degree_map->at({1, 1})->insert({r, r_vertex->get_degree()});
            right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
            right_index_map->at(r)->insert(1, 1);

            right_vertex_mutex_map->insert({r, make_shared<mutex>()});
        }

        auto total_mutex = make_shared<mutex>();
        auto total_pair_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
        total_pair_set->insert({1, 1});

        auto current_pair_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair,equal_pair>>();
        auto next_pair_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair,equal_pair>>();
        current_pair_set->insert({1, 1});

        auto pool = make_shared<thread_pool>(thread_number);
        while (!current_pair_set->empty()) {
           for(const auto&p:*current_pair_set){
               pool->submit_task([=] {
                   auto left_vertex_degree_map = core_left_vertex_degree_map->at(p);
                   auto right_vertex_degree_map = core_right_vertex_degree_map->at(p);

                   auto[i, j] = p;

                   total_mutex->lock();
                   if (!total_pair_set->count({i + 1, j}) && !total_pair_set->count({i, j + 1})) {
                       {
                           total_pair_set->insert({i + 1, j});
                           next_pair_set->insert({i + 1, j});
                           total_pair_set->insert({i, j + 1});
                           next_pair_set->insert({i, j + 1});
                           total_mutex->unlock();

                           auto sub_left_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                                   left_vertex_degree_map);
                           auto sub_right_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                                   right_vertex_degree_map);

                           basic_bipartite_core_decomposition::find_left_core(G, sub_left_vertex_degree_map,
                                                                              sub_right_vertex_degree_map, i + 1, j);
                           if (!sub_left_vertex_degree_map->empty()) {
                               for (const auto&[l, l_degree]: *sub_left_vertex_degree_map) {
                                   left_vertex_mutex_map->at(l)->lock();
                                   left_index_map->at(l)->insert(i + 1, j);
                                   left_vertex_mutex_map->at(l)->unlock();
                               }
                               for (const auto&[r, r_degree]: *sub_right_vertex_degree_map) {
                                   right_vertex_mutex_map->at(r)->lock();
                                   right_index_map->at(r)->insert(j, i + 1);
                                   right_vertex_mutex_map->at(r)->unlock();
                               }

                               core_mutex->lock();
                               core_left_vertex_degree_map->insert({make_pair(i + 1, j), sub_left_vertex_degree_map});
                               core_right_vertex_degree_map->insert({make_pair(i + 1, j), sub_right_vertex_degree_map});
                               core_mutex->unlock();
                           }else
                           {
                               total_mutex->lock();
                               total_pair_set->erase({i + 1, j});
                               next_pair_set->erase({i + 1, j});
                               total_mutex->unlock();
                           }

                       }
                       {
                           find_right_core(G, left_vertex_degree_map, right_vertex_degree_map, i, j + 1);
                           if (!left_vertex_degree_map->empty()) {
                               for (const auto&[l, l_degree]: *left_vertex_degree_map) {
                                   left_vertex_mutex_map->at(l)->lock();
                                   left_index_map->at(l)->insert(i, j + 1);
                                   left_vertex_mutex_map->at(l)->unlock();
                               }
                               for (const auto&[r, r_degree]: *right_vertex_degree_map) {
                                   right_vertex_mutex_map->at(r)->lock();
                                   right_index_map->at(r)->insert(j + 1, i);
                                   right_vertex_mutex_map->at(r)->unlock();
                               }

                               core_mutex->lock();
                               core_left_vertex_degree_map->insert({make_pair(i, j + 1), left_vertex_degree_map});
                               core_right_vertex_degree_map->insert({make_pair(i, j + 1), right_vertex_degree_map});
                               core_mutex->unlock();
                           }else
                           {
                               total_mutex->lock();
                               total_pair_set->erase({i, j + 1});
                               next_pair_set->erase({i, j + 1});
                               total_mutex->unlock();
                           }
                       }
                   } else if (!total_pair_set->count({i + 1, j})) {
                       total_pair_set->insert({i + 1, j});
                       next_pair_set->insert({i + 1, j});
                       total_mutex->unlock();

                       find_left_core(G, left_vertex_degree_map, right_vertex_degree_map, i + 1, j);
                       if (!left_vertex_degree_map->empty()) {
                           for (const auto&[l, l_degree]: *left_vertex_degree_map) {
                               left_vertex_mutex_map->at(l)->lock();
                               left_index_map->at(l)->insert(i + 1, j);
                               left_vertex_mutex_map->at(l)->unlock();
                           }
                           for (const auto&[r, r_degree]: *right_vertex_degree_map) {
                               right_vertex_mutex_map->at(r)->lock();
                               right_index_map->at(r)->insert(j, i + 1);
                               right_vertex_mutex_map->at(r)->unlock();
                           }

                           core_mutex->lock();
                           core_left_vertex_degree_map->insert({{i + 1, j}, left_vertex_degree_map});
                           core_right_vertex_degree_map->insert({{i + 1, j}, right_vertex_degree_map});
                           core_mutex->unlock();
                       }else
                       {
                           total_mutex->lock();
                           total_pair_set->erase({i + 1, j});
                           next_pair_set->erase({i + 1, j});
                           total_mutex->unlock();
                       }
                   } else if (!total_pair_set->count({i, j + 1})) {
                       total_pair_set->insert({i, j + 1});
                       next_pair_set->insert({i, j + 1});
                       total_mutex->unlock();

                       find_right_core(G, left_vertex_degree_map, right_vertex_degree_map, i, j + 1);
                       if (!left_vertex_degree_map->empty()) {
                           for (const auto&[l, l_degree]: *left_vertex_degree_map) {
                               left_vertex_mutex_map->at(l)->lock();
                               left_index_map->at(l)->insert(i, j + 1);
                               left_vertex_mutex_map->at(l)->unlock();
                           }
                           for (const auto&[r, r_degree]: *right_vertex_degree_map) {
                               right_vertex_mutex_map->at(r)->lock();
                               right_index_map->at(r)->insert(j + 1, i);
                               right_vertex_mutex_map->at(r)->unlock();
                           }

                           core_mutex->lock();
                           core_left_vertex_degree_map->insert({make_pair(i, j + 1), left_vertex_degree_map});
                           core_right_vertex_degree_map->insert({make_pair(i, j + 1), right_vertex_degree_map});
                           core_mutex->unlock();
                       }else
                       {
                           total_mutex->lock();
                           total_pair_set->erase({i, j + 1});
                           next_pair_set->erase({i, j + 1});
                           total_mutex->unlock();
                       }
                   }else
                   {
                       total_mutex->unlock();
                   }
              });
           }
           current_pair_set->clear();
           pool->barrier();
           swap(*current_pair_set,*next_pair_set);
        }

        return total_pair_set;
    }



    void basic_bipartite_core_decomposition::find_left_core(const shared_ptr<abstract_bipartite_graph> &G,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_degree_map,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_degree_map,
                                                             uint32_t i,
                                                             uint32_t j) {
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();
        for (const auto&[l, l_degree]: *left_degree_map) {
            if (l_degree < i) {
                l_set->insert(l);
            }
        }

        remove_unsatisfied_vertices(G,left_degree_map,right_degree_map,l_set, r_set,i,j);
    }

    void basic_bipartite_core_decomposition::find_right_core(const shared_ptr<abstract_bipartite_graph> &G,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_degree_map,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_degree_map,
                                                              uint32_t i,
                                                              uint32_t j) {
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &[r, r_degree]: *right_degree_map) {
            if (r_degree < j) {
                r_set->insert(r);
            }
        }
        remove_unsatisfied_vertices(G,left_degree_map,right_degree_map,l_set, r_set,i,j);
    }

    void basic_bipartite_core_decomposition::remove_unsatisfied_vertices(const shared_ptr<abstract_bipartite_graph> &G,
                                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &left_degree_map,
                                                                          const shared_ptr<unordered_map<uint32_t, uint32_t>> &right_degree_map,
                                                                          const shared_ptr<unordered_set<uint32_t>> &l_set,
                                                                          const shared_ptr<unordered_set<uint32_t>> &r_set,
                                                                          uint32_t i,
                                                                          uint32_t j) {
        while (!l_set->empty() || !r_set->empty()) {
            while (!l_set->empty()) {
                auto l = *l_set->begin();
                l_set->erase(l);
                left_degree_map->erase(l);

                for (const auto&[r, e]: *G->get_left_vertex(l)->get_edge_map()) {
                    if (right_degree_map->count(r) && right_degree_map->at(r) >= j) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) < j) {
                            r_set->insert(r);
                        }
                    }
                }
            }

            while (!r_set->empty()) {
                auto r = *r_set->begin();
                r_set->erase(r);
                right_degree_map->erase(r);

                for (const auto &[l, e]: *G->get_right_vertex(r)->get_edge_map()) {
                    if (left_degree_map->count(l) && left_degree_map->at(l) >= i) {
                        --left_degree_map->at(l);
                        if (left_degree_map->at(l) < i) {
                            l_set->insert(l);
                        }
                    }
                }
            }
        }
    }
}

