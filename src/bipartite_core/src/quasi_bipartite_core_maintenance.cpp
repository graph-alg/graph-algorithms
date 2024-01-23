
#include "bipartite_core/quasi_bipartite_core_maintenance.h"

namespace scnu{
    uint32_t quasi_bipartite_core_maintenance::compute_left_vertex_core_degree(
            const shared_ptr<abstract_bipartite_graph> &G,
            const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
            uint32_t l,
            uint32_t i,
            uint32_t j)
    {
        auto degree = 0;
        for(const auto&[r,e]:*G->get_left_vertex(l)->get_edge_map()){
            if(right_index_map->at(r)->count(j, i)){
                ++degree;
            }
        }
        return degree;
    }

    uint32_t
    quasi_bipartite_core_maintenance::compute_right_vertex_core_degree(const shared_ptr<abstract_bipartite_graph> &G,
                                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                                       uint32_t r,
                                                                       uint32_t i,
                                                                       uint32_t j)
    {
        auto degree = 0;
        for(const auto&[l,e]:*G->get_right_vertex(r)->get_edge_map()){
            if(left_index_map->at(l)->count(i,j)){
                ++degree;
            }
        }
        return degree;
    }

    bool quasi_bipartite_core_maintenance::left_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                           uint32_t i,
                                                           uint32_t j,
                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_l_map,
                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_r_map) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto left_core_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_core_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t count = 0;
        for(const auto&e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();
            auto max_j = left_index_map->at(l)->get_j(i);
            auto max_i = right_index_map->at(r)->get_i(j);
            if(max_j < j - 1 || max_i < i -1){
                ++count;
            }else{
                if (max_j == j - 1) {
                    if (!left_core_degree_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, i - 1, j);
                        left_core_degree_map->insert({l, l_degree});
                    }

                    if (left_core_degree_map->at(l) >= i) {
                        l_set->insert(l);
                    } else {
                        evicted_l_set->insert(l);
                    }
                }

                if (max_i == i - 1) {
                    r_set->insert(r);
                }
            }
        }

        if(count == edge_set->size()){
            return false;
        }

        while (!l_set->empty() || !r_set->empty()) {
            while (!l_set->empty()) {
                auto l = *l_set->begin();
                l_set->erase(l);

                if(evicted_l_set->count(l)){
                    continue;
                }

                auto affected_r_set = make_shared<unordered_set<uint32_t> >();

                auto degree  = 0;
                for (const auto& [r,e]:*B->get_left_vertex(l)->get_edge_map()) {
                    auto max_i = right_index_map->at(r)->get_i(j);
                    if(max_i < i -1){
                        continue;
                    }
                    if(max_i > i - 1 || candidate_r_map->count(r)){
                        ++degree;
                    }else if(!evicted_r_set->count(r)){
                        ++degree;
                        affected_r_set->insert(r);
                    }
                }

                if (degree >= i) {
                    candidate_l_map->insert({l, degree});
                    r_set->merge(*affected_r_set);
                } else {
                    remove_left_vertex(B,candidate_l_map,candidate_r_map,evicted_l_set, evicted_r_set, l, i, j);
                }
            }

            while (!r_set->empty()) {
                auto r = *r_set->begin();
                r_set->erase(r);

                if(evicted_r_set->count(r)){
                    continue;
                }

                auto affected_l_set = make_shared<unordered_set<uint32_t> >();

                auto degree = 0;
                auto r_vertex = B->get_right_vertex(r);
                for (const auto& [l,e]:*r_vertex->get_edge_map()) {
                    auto max_j = left_index_map->at(l)->get_j(i);
                    if(max_j < j - 1){
                        continue;
                    }
                    if(max_j > j - 1||candidate_l_map->count(l)){
                        ++degree;
                    }else if(!evicted_l_set->count(l)){
                        if (!left_core_degree_map->count(l)) {
                            auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, i - 1, j);
                            left_core_degree_map->insert({l, l_degree});
                        }
                        if (left_core_degree_map->at(l) >= i) {
                            ++degree;
                            affected_l_set->insert(l);
                        } else {
                            evicted_l_set->insert(l);
                        }
                    }
                }

                if (degree >= j) {
                    candidate_r_map->insert({r, degree});
                    l_set->merge(*affected_l_set);
                } else {
                    remove_right_vertex(B,candidate_l_map,candidate_r_map,evicted_l_set, evicted_r_set, r, i, j);
                }
            }
        }
        return true;
    }

    bool quasi_bipartite_core_maintenance::right_candidate_graph(const shared_ptr<abstract_bipartite_graph> &B,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                           const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                           uint32_t i,
                                                           uint32_t j,
                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_l_map,
                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> & candidate_r_map) {
        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        auto left_core_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_core_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t count = 0;
        for(const auto&e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();
            auto max_i = right_index_map->at(r)->get_i(j);
            auto max_j = left_index_map->at(l)->get_j(i);
            if(max_i < i -1 || max_j < j -1){
                ++count;
            }else {
                if (max_i == i - 1) {
                    l_set->insert(l);
                }

                if (max_j == j - 1) {
                    if (!right_core_degree_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, i, j - 1);
                        right_core_degree_map->insert({r, r_degree});
                    }

                    if (right_core_degree_map->at(r) >= j) {
                        r_set->insert(r);
                    } else {
                        evicted_r_set->insert(r);
                    }
                }
            }
        }

        if(count == edge_set->size()){
            return false;
        }

        while (!l_set->empty() || !r_set->empty()) {
            while (!l_set->empty()) {
                auto l = *l_set->begin();
                l_set->erase(l);

                if(evicted_l_set->count(l)){
                    continue;
                }

                auto affected_r_set = make_shared<unordered_set<uint32_t> >();

                auto degree  = 0;
                for (const auto& [r,e]:*B->get_left_vertex(l)->get_edge_map()) {
                    auto max_i = right_index_map->at(r)->get_i(j);
                    if(max_i < i -1){
                        continue;
                    }
                    if(max_i > i - 1 || candidate_r_map->count(r)){
                        ++degree;
                    }else if(!evicted_r_set->count(r)){
                        if (!right_core_degree_map->count(r)) {
                            auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, i, j - 1);
                            right_core_degree_map->insert({r, r_degree});
                        }
                        if (right_core_degree_map->at(r) >= j) {
                            ++degree;
                            affected_r_set->insert(r);
                        } else {
                            evicted_r_set->insert(r);
                        }
                    }
                }

                if (degree >= i) {
                    candidate_l_map->insert({l, degree});
                    r_set->merge(*affected_r_set);
                } else {
                    remove_left_vertex(B,candidate_l_map,candidate_r_map,evicted_l_set, evicted_r_set, l, i, j);
                }
            }

            while (!r_set->empty()) {
                auto r = *r_set->begin();
                r_set->erase(r);

                if(evicted_r_set->count(r)){
                    continue;
                }

                auto affected_l_set = make_shared<unordered_set<uint32_t> >();

                auto degree = 0;
                auto r_vertex = B->get_right_vertex(r);
                for (const auto& [l,e]:*r_vertex->get_edge_map()) {
                    auto max_j = left_index_map->at(l)->get_j(i);
                    if(max_j < j -1){
                        continue;
                    }
                    if(max_j > j - 1 || candidate_l_map->count(l)) {
                        ++degree;
                    }else if(!evicted_l_set->count(l)){
                        ++degree;
                        affected_l_set->insert(l);
                    }
                }

                if (degree >= j) {
                    candidate_r_map->insert({r, degree});
                    l_set->merge(*affected_l_set);
                } else {
                    remove_right_vertex(B,candidate_l_map,candidate_r_map,evicted_l_set, evicted_r_set, r, i, j);
                }
            }
        }
        return true;
    }

    void quasi_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()){
            new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
        }
    }

    void quasi_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                const shared_ptr<thread_pool>& pool) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, shared_ptr<bipartite_core_left_store_index>()});
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            new_right_index_map->insert({r, shared_ptr<bipartite_core_right_store_index>()});
        }

        for(const auto &[l, l_index]:*new_left_index_map){
            pool->submit_task([=]{
                new_left_index_map->at(l) =  make_shared<bipartite_core_left_store_index>();
            });
        }

        for(const auto &[r, r_index]:*new_right_index_map){
            pool->submit_task([=] {
                new_right_index_map->at(r) = make_shared<bipartite_core_right_store_index>();
            });
        }
        pool->barrier();
    }


    void quasi_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map) {
        {
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
        }

        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto quasi_l_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        auto quasi_r_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        {
            quasi_l_map->insert({{2, 2}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            quasi_r_map->insert({{2, 2}, make_shared<unordered_map<uint32_t, uint32_t>>()});

            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                quasi_l_map->at({2, 2})->insert({l, B->get_left_vertex(l)->get_degree()});

                auto r = e->get_right_vertex_id();
                quasi_r_map->at({2, 2})->insert({r, B->get_right_vertex(r)->get_degree()});
            }
            auto sub_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(edge_set);
            find_quasi_core(B, sub_edge_set, quasi_l_map->at({2,2}), quasi_r_map->at({2,2}), 2, 2);
            quasi_edge_map->insert({{2, 2}, sub_edge_set});
        }

        while(!quasi_edge_map->empty())
        {
            auto flag_set = make_shared<unordered_set<pair<uint32_t,uint32_t>, hash_pair, equal_pair>>();
            for(const auto&p:*quasi_edge_map)
            {
                auto [i, j] = p.first;

                auto sub_edge_set = p.second;

                auto candidate_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                auto candidate_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                if(i == 2){
                    auto flag = right_candidate_graph(B, sub_edge_set, left_index_map, right_index_map, i, j,
                                                candidate_l_map, candidate_r_map);

                    if (flag) {
                        for (const auto &[l, l_degree]: *candidate_l_map) {
                            new_left_index_map->at(l)->insert(i, j);
                        }
                        for (const auto &[r, r_degree]: *candidate_r_map) {
                            new_right_index_map->at(r)->insert(j, i);
                        }
                    }else{
                        flag_set->insert(p.first);
                    }
                }
                else{
                    auto flag = left_candidate_graph(B, sub_edge_set, left_index_map, right_index_map, i, j,
                                                candidate_l_map, candidate_r_map);

                    if (flag) {
                        for (const auto &[l, l_degree]: *candidate_l_map) {
                            new_left_index_map->at(l)->insert(i, j);
                        }
                        for (const auto &[r, r_degree]: *candidate_r_map) {
                            new_right_index_map->at(r)->insert(j, i);
                        }
                    }else{
                        flag_set->insert(p.first);
                    }
                }
            }
            while(!quasi_edge_map->empty()){
                auto p = *quasi_edge_map->begin();
                auto [i, j] = p.first;
                auto sub_edge_set = p.second;
                quasi_edge_map->erase(p.first);

                get_quasi_cores(B, sub_edge_set, quasi_l_map,
                                    quasi_r_map, i, j, next_quasi_edge_map);
            }
            swap(*quasi_edge_map, *next_quasi_edge_map);
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
            for(const auto &p:*flag_set){
                auto [i, j] = p;
                quasi_edge_map->erase({i + 1, j});
                quasi_edge_map->erase({i, j + 1});
            }
        }
    }

    void quasi_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                  const shared_ptr<thread_pool>& pool) {
        {
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
        }

        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto quasi_l_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        auto quasi_r_map = make_shared<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>>();
        {
            quasi_l_map->insert({{2, 2}, make_shared<unordered_map<uint32_t, uint32_t>>()});
            quasi_r_map->insert({{2, 2}, make_shared<unordered_map<uint32_t, uint32_t>>()});

            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                quasi_l_map->at({2, 2})->insert({l, B->get_left_vertex(l)->get_degree()});

                auto r = e->get_right_vertex_id();
                quasi_r_map->at({2, 2})->insert({r, B->get_right_vertex(r)->get_degree()});
            }
            auto sub_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(edge_set);
            find_quasi_core(B, sub_edge_set, quasi_l_map->at({2,2}), quasi_r_map->at({2,2}), 2, 2);
            quasi_edge_map->insert({{2, 2}, sub_edge_set});
        }

        auto global_left_index_mutex = make_shared<mutex>();
        auto global_right_index_mutex = make_shared<mutex>();

        while(!quasi_edge_map->empty())
        {
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            for(const auto&p:*quasi_edge_map)
            {
                pool->submit_task([=]{
                    auto [i, j] = p.first;
                    auto sub_edge_set = p.second;

                    auto candidate_l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    auto candidate_r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

                    if(i == 2){
                        auto flag = right_candidate_graph(B, sub_edge_set, left_index_map, right_index_map, i, j,
                                                    candidate_l_map, candidate_r_map);

                        if (flag) {
                            global_left_index_mutex->lock();
                            for (const auto &[l, l_degree]: *candidate_l_map) {
                                new_left_index_map->at(l)->insert(i, j);
                            }
                            global_left_index_mutex->unlock();

                            global_right_index_mutex->lock();
                            for (const auto &[r, r_degree]: *candidate_r_map) {
                                new_right_index_map->at(r)->insert(j, i);
                            }
                            global_right_index_mutex->unlock();

                        } else {
                            global_left_index_mutex->lock();
                            flag_set->insert(p.first);
                            global_left_index_mutex->unlock();
                        }
                    }else{
                        auto flag = left_candidate_graph(B, sub_edge_set, left_index_map, right_index_map, i, j,
                                                          candidate_l_map, candidate_r_map);
                        if (flag) {
                            global_left_index_mutex->lock();
                            for (const auto &[l, l_degree]: *candidate_l_map) {
                                new_left_index_map->at(l)->insert(i, j);
                            }
                            global_left_index_mutex->unlock();

                            global_right_index_mutex->lock();
                            for (const auto &[r, r_degree]: *candidate_r_map) {
                                new_right_index_map->at(r)->insert(j, i);
                            }
                            global_right_index_mutex->unlock();
                        } else {
                            global_left_index_mutex->lock();
                            flag_set->insert(p.first);
                            global_left_index_mutex->unlock();
                        }
                    }
                });
            }
            while(!quasi_edge_map->empty()){
                auto p = *quasi_edge_map->begin();
                auto [i, j] = p.first;
                auto sub_edge_set = p.second;
                quasi_edge_map->erase(p.first);

                get_quasi_cores(B, sub_edge_set, quasi_l_map,
                                quasi_r_map, i, j, next_quasi_edge_map);
            }
            swap(*quasi_edge_map, *next_quasi_edge_map);
            pool->barrier();

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
            for(const auto &p:*flag_set){
                auto [i, j] = p;
                quasi_edge_map->erase({i + 1, j});
                quasi_edge_map->erase({i, j + 1});
            }
        }
    }

    void quasi_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map)
    {
        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto quasi_l_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_map<uint32_t,uint32_t>>,hash_pair,equal_pair>>();
        auto quasi_r_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_map<uint32_t,uint32_t>>,hash_pair,equal_pair>>();
        {

            quasi_l_map->insert({{2, 2}, make_shared<unordered_map<uint32_t,uint32_t>>()});
            quasi_r_map->insert({{2, 2}, make_shared<unordered_map<uint32_t,uint32_t>>()});

            for(const auto&e:*edge_set)
            {
                auto l = e->get_left_vertex_id();
                quasi_l_map->at({2, 2})->insert({l, B->get_left_vertex(l)->get_degree()});

                auto r = e->get_right_vertex_id();
                quasi_r_map->at({2, 2})->insert({r, B->get_right_vertex(r)->get_degree()});
            }

            auto sub_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(edge_set);
            find_quasi_core(B, sub_edge_set, quasi_l_map->at({2,2}), quasi_r_map->at({2,2}), 2, 2);
            quasi_edge_map->insert({{2, 2}, sub_edge_set});
        }

        auto isolated_l_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_r_set = make_shared<unordered_set<uint32_t>>();
        B->remove_edge_collection(edge_set, isolated_l_set, isolated_r_set);

        while (!quasi_edge_map->empty())
        {
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            for(const auto &p:*quasi_edge_map){
                auto [i, j] = p.first;

                auto sub_edge_set = p.second;
                auto sub_removed_l_set = make_shared<unordered_set<uint32_t>>();
                auto sub_removed_r_set = make_shared<unordered_set<uint32_t>>();

                auto flag = update_single_core(B, sub_edge_set, left_index_map, right_index_map, sub_removed_l_set,
                                               sub_removed_r_set, i, j);

                if (flag) {
                    for (const auto &l: *sub_removed_l_set) {
                        new_left_index_map->at(l)->remove(i, j);
                    }
                    for (const auto &r: *sub_removed_r_set) {
                        new_right_index_map->at(r)->remove(j, i);
                    }
                }else{
                    flag_set->insert(p.first);
                }
            }
            B->insert_edge_collection(edge_set);
            while(!quasi_edge_map->empty()){
                auto p = *quasi_edge_map->begin();
                auto [i, j] = p.first;
                auto sub_edge_set = p.second;
                quasi_edge_map->erase(p.first);

                if(!flag_set->count(p.first)){
                    get_quasi_cores(B, sub_edge_set, quasi_l_map, quasi_r_map, i, j, next_quasi_edge_map);
                }
            }
            B->remove_edge_collection(edge_set);
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }

        for(const auto&[l, l_index]:*new_left_index_map){
           for(const auto &[i, j]:*l_index->get_index_map()){
               left_index_map->at(l)->remove(i, j - 1);
               if(left_index_map->at(l)->get_j(i) == 0){
                   left_index_map->at(l)->remove(i);
               }
           }
           l_index->clear();
        }
        for(const auto&[r, r_index]:*new_right_index_map){
            for(const auto &[j, i]:*r_index->get_index_map()){
                right_index_map->at(r)->remove(j, i - 1);
                if(right_index_map->at(r)->get_i(j) == 0){
                    right_index_map->at(r)->remove(j);
                }
            }
            r_index->clear();
        }

        update_trivial_bipartite_cores(B, edge_set, left_index_map, right_index_map);

        for(const auto&l:*isolated_l_set)
        {
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }

        for(const auto&r:*isolated_r_set){
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }
    }

    void quasi_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                  const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& new_left_index_map,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& new_right_index_map,
                                                  const shared_ptr<thread_pool>& pool)
    {
        auto quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto next_quasi_edge_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>>();
        auto quasi_l_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_map<uint32_t,uint32_t>>,hash_pair,equal_pair>>();
        auto quasi_r_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,shared_ptr<unordered_map<uint32_t,uint32_t>>,hash_pair,equal_pair>>();
        {

            quasi_l_map->insert({{2, 2}, make_shared<unordered_map<uint32_t,uint32_t>>()});
            quasi_r_map->insert({{2, 2}, make_shared<unordered_map<uint32_t,uint32_t>>()});

            for(const auto&e:*edge_set)
            {
                auto l = e->get_left_vertex_id();
                quasi_l_map->at({2, 2})->insert({l, B->get_left_vertex(l)->get_degree()});

                auto r = e->get_right_vertex_id();
                quasi_r_map->at({2, 2})->insert({r, B->get_right_vertex(r)->get_degree()});
            }

            auto sub_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(edge_set);
            find_quasi_core(B, sub_edge_set, quasi_l_map->at({2,2}), quasi_r_map->at({2,2}), 2, 2);
            quasi_edge_map->insert({{2, 2}, sub_edge_set});
        }

        auto isolated_l_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_r_set = make_shared<unordered_set<uint32_t>>();
        B->remove_edge_collection(edge_set, isolated_l_set, isolated_r_set);

        auto global_left_index_mutex = make_shared<mutex>();
        auto global_right_index_mutex = make_shared<mutex>();
        while (!quasi_edge_map->empty())
        {
            auto flag_set = make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair, equal_pair>>();
            for(const auto &p:*quasi_edge_map){
                pool->submit_task([=]{
                    auto [i, j] = p.first;
                    auto sub_edge_set = p.second;

                    auto sub_removed_l_set = make_shared<unordered_set<uint32_t>>();
                    auto sub_removed_r_set = make_shared<unordered_set<uint32_t>>();

                    auto flag = update_single_core(B, sub_edge_set, left_index_map, right_index_map, sub_removed_l_set,
                                                   sub_removed_r_set, i, j);

                    if (flag) {
                        global_left_index_mutex->lock();
                        for (const auto &l: *sub_removed_l_set) {
                            new_left_index_map->at(l)->remove(i, j);
                        }
                        global_left_index_mutex->unlock();

                        global_right_index_mutex->lock();
                        for (const auto &r: *sub_removed_r_set) {
                            new_right_index_map->at(r)->remove(j, i);
                        }
                        global_right_index_mutex->unlock();
                    }else{
                        global_left_index_mutex->lock();
                        flag_set->insert(p.first);
                        global_left_index_mutex->unlock();
                    }
                });
            }
            pool->barrier();
            B->insert_edge_collection(edge_set);
            while(!quasi_edge_map->empty()){
                auto p = *quasi_edge_map->begin();
                auto [i, j] = p.first;
                auto sub_edge_set = p.second;
                quasi_edge_map->erase(p.first);

                if(!flag_set->count(p.first)){
                    get_quasi_cores(B, sub_edge_set, quasi_l_map, quasi_r_map, i, j, next_quasi_edge_map);
                }
            }
            B->remove_edge_collection(edge_set);
            swap(*quasi_edge_map, *next_quasi_edge_map);
        }

        for(const auto&[l, l_index]:*new_left_index_map){
            for(const auto &[i, j]:*l_index->get_index_map()){
                left_index_map->at(l)->remove(i, j - 1);
                if(left_index_map->at(l)->get_j(i) == 0){
                    left_index_map->at(l)->remove(i);
                }
            }
            l_index->clear();
        }
        for(const auto&[r, r_index]:*new_right_index_map){
            for(const auto &[j, i]:*r_index->get_index_map()){
                right_index_map->at(r)->remove(j, i - 1);
                if(right_index_map->at(r)->get_i(j) == 0){
                    right_index_map->at(r)->remove(j);
                }
            }
            r_index->clear();
        }

        update_trivial_bipartite_cores(B, edge_set, left_index_map, right_index_map, pool);

        for(const auto&l:*isolated_l_set)
        {
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }

        for(const auto&r:*isolated_r_set){
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }
    }

    void quasi_bipartite_core_maintenance::find_quasi_core(const shared_ptr<abstract_bipartite_graph> &G,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                                           const shared_ptr<unordered_map<uint32_t,uint32_t>>& l_map,
                                                           const shared_ptr<unordered_map<uint32_t,uint32_t>>& r_map,
                                                           uint32_t i,
                                                           uint32_t j) {
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();
        //find unsatisfied left vertices
        for (const auto&[l,degree]: *l_map) {
            if (degree < i) {
                l_set->insert(l);
            }
        }
        //find unsatisfied right vertices
        for (const auto& [r,degree]:*r_map) {
            if (degree < j) {
                r_set->insert(r);
            }
        }

        //delete unsatisfied vertices and edges
        removal_unsatisfied_vertices(G,edge_set,l_map,r_map, l_set, r_set,i, j);
    }



    void quasi_bipartite_core_maintenance::get_quasi_cores(const shared_ptr<abstract_bipartite_graph> &B,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                                           const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_core_left_vertex_degree_map,
                                                           const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_map<uint32_t, uint32_t>>, hash_pair, equal_pair>> &quasi_core_right_vertex_degree_map,
                                                           uint64_t i,
                                                           uint64_t j,
                                                           const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>, hash_pair, equal_pair>> &next_quasi_edge_map) {
        auto left_vertex_degree_map = quasi_core_left_vertex_degree_map->at({i, j});
        auto right_vertex_degree_map = quasi_core_right_vertex_degree_map->at({i, j});
        quasi_core_left_vertex_degree_map->erase({i, j});
        quasi_core_right_vertex_degree_map->erase({i, j});

        if (i == 2) {
            auto sub_left_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                    left_vertex_degree_map);
            auto sub_right_vertex_degree_map = container_copy::to_unordered_map<uint32_t, uint32_t>(
                    right_vertex_degree_map);

            auto sub_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    edge_set);
            find_quasi_core(B, sub_edge_set, sub_left_vertex_degree_map,
                            sub_right_vertex_degree_map, i, j + 1);
            if (!sub_edge_set->empty()) {
                quasi_core_left_vertex_degree_map->insert(
                        {make_pair(i, j + 1), sub_left_vertex_degree_map});
                quasi_core_right_vertex_degree_map->insert(
                        {make_pair(i, j + 1), sub_right_vertex_degree_map});
                next_quasi_edge_map->insert({{i, j + 1}, sub_edge_set});
            }
        }

        {
            auto sub_edge_set = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(
                    edge_set);

            find_quasi_core(B, sub_edge_set, left_vertex_degree_map,
                            right_vertex_degree_map, i + 1, j);

            if (!sub_edge_set->empty()) {
                quasi_core_left_vertex_degree_map->insert({make_pair(i + 1, j), left_vertex_degree_map});
                quasi_core_right_vertex_degree_map->insert({make_pair(i + 1, j), right_vertex_degree_map});

                next_quasi_edge_map->insert({{i + 1, j}, sub_edge_set});
            }
        }
    }

    void quasi_bipartite_core_maintenance::removal_unsatisfied_vertices(const shared_ptr<abstract_bipartite_graph> &B,
                                                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &l_map,
                                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &r_map,
                                                                        const shared_ptr<unordered_set<uint32_t>> &l_set,
                                                                        const shared_ptr<unordered_set<uint32_t>> &r_set,
                                                                        uint32_t i,
                                                                        uint32_t j) {
        while (!l_set->empty() || !r_set->empty()) {
            while (!l_set->empty()) {
                auto l = *l_set->begin();
                l_set->erase(l);

                auto l_vertex = B->get_left_vertex(l);
                for (const auto&[r,e]:*l_vertex->get_edge_map()) {
                    edge_set->erase(e);
                    if(r_map->count(r)){
                        --r_map->at(r);
                        if (r_map->at(r)< j) {
                            r_map->erase(r);
                            r_set->insert(r);
                        }
                    }
                }
                l_map->erase(l);
            }

            while (!r_set->empty()) {
                auto r = *r_set->begin();
                r_set->erase(r);

                auto r_vertex = B->get_right_vertex(r);
                for (const auto& [l,e]:*r_vertex->get_edge_map()) {
                    edge_set->erase(e);
                    if(l_map->count(l)){
                        --l_map->at(l);
                        if (l_map->at(l) < i) {
                            l_map->erase(l);
                            l_set->insert(l);
                        }
                    }
                }
                r_map->erase(r);
            }
        }
    }


    void quasi_bipartite_core_maintenance::remove_left_vertex(const shared_ptr<abstract_bipartite_graph>& B,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_l_map,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_r_map,
                                                              const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                              const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                              const uint32_t l,
                                                              uint32_t i,
                                                              uint32_t j) {
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        visited_set->insert(l);

        while(!visited_set->empty()){
            auto l1 = *visited_set->begin();
            visited_set->erase(l1);
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
                                    visited_set->insert(l2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void quasi_bipartite_core_maintenance::remove_right_vertex(const shared_ptr<abstract_bipartite_graph>& B,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_l_map,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>>& candidate_r_map,
                                                              const shared_ptr<unordered_set<uint32_t>> &evicted_l_set,
                                                              const shared_ptr<unordered_set<uint32_t>> &evicted_r_set,
                                                              const uint32_t r,
                                                              uint32_t i,
                                                              uint32_t j)
    {
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        visited_set->insert(r);

        while(!visited_set->empty()){
            auto r1 = *visited_set->begin();
            visited_set->erase(r1);
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
                                    visited_set->insert(r2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    uint32_t quasi_bipartite_core_maintenance::update_single_core(const shared_ptr<abstract_bipartite_graph>& B,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                              const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>>& left_index_map,
                                                              const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>>& right_index_map,
                                                              const shared_ptr<unordered_set<uint32_t>>& removed_l_set,
                                                              const shared_ptr<unordered_set<uint32_t>>& removed_r_set,
                                                              uint32_t i,
                                                              uint32_t j) {
        auto l_set = make_shared<unordered_set<uint32_t>>();
        auto r_set = make_shared<unordered_set<uint32_t>>();

        auto l_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto r_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        uint32_t count = 0;
        for (const auto &e: *edge_set) {
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if (left_index_map->at(l)->count(i, j) && right_index_map->at(r)->count(j, i)) {
                if(B->get_left_vertex(l)){
                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, i, j);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) < i) {
                        l_set->insert(l);
                    }
                }
                if(B->get_right_vertex(r)){
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, i, j);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) < j) {
                        r_set->insert(r);
                    }
                }

            }
            else
            {
                ++count;
            }
        }

        if(count == edge_set->size()){
            return false;
        }
        while (!l_set->empty() || !r_set->empty()) {
            while (!l_set->empty()) {
                auto l = *l_set->begin();
                l_set->erase(l);
                removed_l_set->insert(l);

                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (removed_r_set->count(r) || !right_index_map->at(r)->count(j, i)) {
                        continue;
                    }
                    if (!r_map->count(r)) {
                        auto r_degree = compute_right_vertex_core_degree(B, left_index_map, r, i, j);
                        r_map->insert({r, r_degree});
                    }

                    if (r_map->at(r) >= j) {
                        --r_map->at(r);
                        if (r_map->at(r) < j) {
                            r_set->insert(r);
                        }
                    }
                }
            }

            while (!r_set->empty()) {
                auto r = *r_set->begin();
                r_set->erase(r);

                removed_r_set->insert(r);
                for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                    if (removed_l_set->count(l) || !left_index_map->at(l)->count(i, j)) {
                        continue;
                    }

                    if (!l_map->count(l)) {
                        auto l_degree = compute_left_vertex_core_degree(B, right_index_map, l, i, j);
                        l_map->insert({l, l_degree});
                    }

                    if (l_map->at(l) >= i) {
                        --l_map->at(l);
                        if (l_map->at(l) < i) {
                            l_set->insert(l);
                        }
                    }
                }
            }
        }
        return true;
    }

    void quasi_bipartite_core_maintenance::update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
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
            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (B->get_right_vertex(r)->get_degree() > max_degree) {
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
            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if (B->get_left_vertex(l)->get_degree() > max_degree) {
                    max_degree = B->get_left_vertex(l)->get_degree();
                }
            }
            right_index_map->at(r)->remove(1, max_degree);
        }
    }

    void quasi_bipartite_core_maintenance::update_trivial_bipartite_cores(const shared_ptr<abstract_bipartite_graph> &B,
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

}

