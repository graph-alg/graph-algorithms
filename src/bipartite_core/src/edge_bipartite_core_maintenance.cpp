
#include "bipartite_core/edge_bipartite_core_maintenance.h"

namespace scnu{
    void edge_bipartite_core_maintenance::add_left_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                              const shared_ptr<unordered_set<uint32_t>> &Tl,
                                                              const shared_ptr<unordered_set<uint32_t>> &Tr,
                                                              const shared_ptr<unordered_set<uint32_t>> &Cl,
                                                              const shared_ptr<unordered_set<uint32_t>> &Cr,
                                                              const shared_ptr<unordered_set<uint32_t>> &Sl,
                                                              const shared_ptr<unordered_set<uint32_t>> &Sr,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_l,
                                                              const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_r,
                                                              uint32_t l,
                                                              uint32_t alpha,
                                                              uint32_t beta)
    {
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        visited_set->insert(l);

        while(!visited_set->empty()){
            auto l1 = *visited_set->begin();
            visited_set->erase(l1);

            Cl->insert(l1);
            for(const auto &[r1,e1]:* B->get_left_vertex(l1)->get_edge_map()){
                if(!Tr->count(r1) && right_index_map->at(r1)->get_i(beta) == alpha){
                    Sr->insert(r1);
                }else if(Tr->count(r1) && !Cr->count(r1) && sup_r->at(r1) >= beta){
                    --sup_r->at(r1);
                    if(sup_r->at(r1) < beta)
                    {
                        Cr->insert(r1);
                        for(const auto &[l2,e2]:*B->get_right_vertex(r1)->get_edge_map()){
                            if(!Tl->count(l2) && left_index_map->at(l2)->get_j(alpha) == beta){
                                Sl->insert(l2);
                            }else if(Tl->count(l2) && !Cl->count(l2) && sup_l->at(l2) >= alpha){
                                --sup_l->at(l2);
                                if(sup_l->at(l2) < alpha){
                                    visited_set->insert(l2);
                                }
                            }
                        }
                    }
                }
            }
        }


    }

    void edge_bipartite_core_maintenance::add_right_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                               const shared_ptr<unordered_set<uint32_t>> &Tl,
                                                               const shared_ptr<unordered_set<uint32_t>> &Tr,
                                                               const shared_ptr<unordered_set<uint32_t>> &Cl,
                                                               const shared_ptr<unordered_set<uint32_t>> &Cr,
                                                               const shared_ptr<unordered_set<uint32_t>> &Sl,
                                                               const shared_ptr<unordered_set<uint32_t>> &Sr,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_l,
                                                               const shared_ptr<unordered_map<uint32_t, uint32_t>> &sup_r,
                                                               uint32_t r,
                                                               uint32_t alpha,
                                                               uint32_t beta)
    {
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        visited_set->insert(r);

        while(!visited_set->empty()){
            auto r1 = *visited_set->begin();
            visited_set->erase(r1);
            Cr->insert(r1);

            for(const auto &[l1,e1]:*B->get_right_vertex(r1)->get_edge_map()){
                if(!Tl->count(l1) && left_index_map->at(l1)->get_j(alpha) == beta){
                    Sl->insert(l1);
                }else if(Tl->count(l1) && !Cl->count(l1) && sup_l->at(l1) >= alpha){
                    --sup_l->at(l1);
                    if(sup_l->at(l1) < alpha){
                        Cl->insert(l1);
                        for(const auto &[r2,e2]:* B->get_left_vertex(l1)->get_edge_map()) {
                            if (!Tr->count(r2) && right_index_map->at(r2)->get_i(beta) == alpha) {
                                Sr->insert(r2);
                            } else if (Tr->count(r2) && !Cr->count(r2) && sup_r->at(r2) >= beta) {
                                --sup_r->at(r2);
                                if (sup_r->at(r2) < beta) {
                                    visited_set->insert(r2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void edge_bipartite_core_maintenance::batch_insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                       const shared_ptr<uint32_t>& delta){
        B->insert_edge_collection(edge_set);
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();
            if(!new_left_index_map->count(l))
            {
                new_left_index_map->insert({l,make_shared<bipartite_core_left_store_index>()});
                new_left_index_map->at(l)->insert(1, 1);

                left_index_map->insert({l,make_shared<bipartite_core_left_store_index>()});
                left_index_map->at(l)->insert(1, 1);
            }
            if(!new_right_index_map->count(r)){
                new_right_index_map->insert({r,make_shared<bipartite_core_right_store_index>()});
                new_right_index_map->at(r)->insert(1, 1);

                right_index_map->insert({r,make_shared<bipartite_core_right_store_index>()});
                right_index_map->at(r)->insert(1, 1);
            }
        }

        *delta = find_core(B);

        for(uint32_t alpha = 1; alpha <= *delta;++alpha){
            uint32_t phi_alpha = UINT32_MAX;
            for(const auto &e:*edge_set){
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                auto min_value = min(get_maximal_b_alpha(B, right_index_map,l, alpha), right_index_map->at(r)->get_maximal_j(alpha));

                if(min_value < phi_alpha){
                    phi_alpha = min_value;
                }
            }

            batch_update_alpha_phi_alpha_core(B, new_left_index_map, new_right_index_map, alpha, phi_alpha);
        }

        for(uint32_t beta = 1; beta <= *delta;++beta){
            uint32_t phi_beta = UINT32_MAX;
            for(const auto &e:*edge_set){
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                auto min_value = min(left_index_map->at(l)->get_maximal_i(beta), get_maximal_b_beta(B, left_index_map, r, beta));

                if(min_value < phi_beta){
                    phi_beta = min_value;
                }
            }

            batch_update_phi_beta_beta_core(B, new_left_index_map, new_right_index_map, phi_beta, beta);
        }

        for(const auto &[l,l_index]:*new_left_index_map){
            for(const auto[i, j]:*l_index->get_index_map()){
                left_index_map->at(l)->insert(i,j);
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

    void edge_bipartite_core_maintenance::batch_insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                       const shared_ptr<uint32_t> &delta,
                                                       const shared_ptr<thread_pool> &pool){

        B->insert_edge_collection(edge_set);
        for(const auto &e:*edge_set){
            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();
            if(!new_left_index_map->count(l))
            {
                new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                left_mutex_map->insert({l, make_shared<mutex>()});

                left_index_map->insert({l,make_shared<bipartite_core_left_store_index>()});
                for(uint32_t i = 1; i <= B->get_left_vertex(l)->get_degree();++i){
                    left_index_map->at(l)->insert(i, 1);
                }

            }
            if(!new_right_index_map->count(r)){
                new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                right_mutex_map->insert({r, make_shared<mutex>()});

                right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                for (uint32_t j = 1; j <= B->get_right_vertex(r)->get_degree(); ++j) {
                    right_index_map->at(r)->insert(j, 1);
                }
            }
        }

        *delta = find_core(B);

        auto alpha_vector = make_shared<vector<uint32_t>>(*delta + 1, 0);
        for (uint32_t alpha = 1; alpha <= *delta; ++alpha) {
            pool->submit_task([=] {
                uint32_t phi_alpha = UINT32_MAX;
                for (const auto &e: *edge_set) {
                    auto l = e->get_left_vertex_id();
                    auto r = e->get_right_vertex_id();

                    auto min_value = std::min(left_index_map->at(l)->get_j(alpha),
                                              right_index_map->at(r)->get_maximal_j(alpha));
                    if (min_value < phi_alpha) {
                        phi_alpha = min_value;
                    }
                }
                alpha_vector->at(alpha) = phi_alpha;
            });
        }

        auto beta_vector = make_shared<vector<uint32_t>>(*delta + 1, 0);
        for (uint32_t beta = 1; beta <= *delta; ++beta) {
            pool->submit_task([=] {
                uint32_t phi_beta = UINT32_MAX;
                for (const auto &e: *edge_set) {
                    auto l = e->get_left_vertex_id();
                    auto r = e->get_right_vertex_id();

                    auto min_value = std::min(left_index_map->at(l)->get_maximal_i(beta),
                                              right_index_map->at(r)->get_i(beta));
                    if (min_value < phi_beta) {
                        phi_beta = min_value;
                    }
                }
                beta_vector->at(beta) = phi_beta;
            });
        }
        pool->barrier();

        for (uint32_t alpha = 1; alpha <= *delta; ++alpha) {
            pool->submit_task([=] {
                auto phi_alpha = alpha_vector->at(alpha);
                batch_update_alpha_phi_alpha_core(B, left_mutex_map, right_mutex_map, new_left_index_map,
                                                  new_right_index_map, alpha, phi_alpha);
            });
        }

        for (uint32_t beta = 1; beta <= *delta; ++beta) {
            pool->submit_task([=]{
                auto phi_beta = beta_vector->at(beta);
                batch_update_phi_beta_beta_core(B, left_mutex_map, right_mutex_map, new_left_index_map,
                                                new_right_index_map, phi_beta, beta);
            });
        }
        pool->barrier();


        for(const auto &[l,l_index]:*new_left_index_map){
            for(const auto[i, j]:*l_index->get_index_map()){
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

    void edge_bipartite_core_maintenance::batch_remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                       const shared_ptr<uint32_t>& delta){
        auto isolated_left_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_right_vertex_set = make_shared<unordered_set<uint32_t>>();
        B->remove_edge_collection(edge_set, isolated_left_vertex_set, isolated_right_vertex_set);

        auto alpha_vector = make_shared<vector<uint32_t>>(*delta + 1,0);
        for(uint32_t alpha = 1;alpha <= *delta; ++alpha){
            uint32_t tau_alpha = 0;
            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                auto max_value = std::min(left_index_map->at(l)->get_j(alpha),
                                          right_index_map->at(r)->get_maximal_j(alpha));
                if (max_value > tau_alpha) {
                    tau_alpha = max_value;
                }
            }
            alpha_vector->at(alpha) = tau_alpha;
        }

        auto beta_vector = make_shared<vector<uint32_t>>(*delta + 1, 0);
        for(uint32_t beta = 1; beta <= *delta; ++beta){
            uint32_t tau_beta = 0;
            for (const auto &e: *edge_set) {
                auto l = e->get_left_vertex_id();
                auto r = e->get_right_vertex_id();

                auto max_value = std::min(left_index_map->at(l)->get_maximal_i(beta),
                                          right_index_map->at(r)->get_i(beta));
                if (max_value > tau_beta) {
                    tau_beta = max_value;
                }
            }
            beta_vector->at(beta) = tau_beta;
        }

        for(uint32_t alpha = 1;alpha <= *delta; ++alpha){
            if(alpha_vector->at(alpha) == 0){
                continue;
            }
            batch_update_alpha_tau_alpha_core(B, new_left_index_map, new_right_index_map, alpha,alpha_vector->at(alpha));
        }

        for(uint32_t beta = 1; beta <= *delta; ++beta){
            if(beta_vector->at(beta) == 0){
                continue;
            }
            batch_update_tau_beta_beta_core(B,new_left_index_map,new_right_index_map, beta_vector->at(beta), beta);
        }

        for(const auto &[l, l_index]:*left_index_map){
            auto index_map = l_index->get_index_map();
            for(const auto &[i,j]:*index_map){
                bool flag  = false;
                for(uint32_t alpha = 1;alpha <= *delta; ++alpha){
                    if(alpha_vector->at(alpha) == 0){
                        continue;
                    }
                    if(i <= alpha && j <= alpha_vector->at(alpha)){
                        flag = true;
                        break;
                    }
                }
                if(!flag){
                    for(uint32_t beta = 1; beta <= *delta; ++beta){
                        if(beta_vector->at(beta) == 0){
                            continue;
                        }
                        if(i <= beta_vector->at(beta) && j <=beta){
                            flag = true;
                            break;
                        }
                    }
                }
                if(!flag){
                    new_left_index_map->at(l)->insert(i,j);
                }
            }
            l_index->clear();
        }

        for(const auto &[r, r_index]:*right_index_map){
            auto index_map = r_index->get_index_map();
            for(const auto &[j,i]:*index_map){
                bool flag  = false;
                for(uint32_t alpha = 1;alpha <= *delta ; ++alpha){
                    if(alpha_vector->at(alpha) == 0){
                        continue;
                    }
                    if(i <= alpha && j <= alpha_vector->at(alpha)){
                        flag = true;
                        break;
                    }
                }
                if(!flag){
                    for(uint32_t beta = 1; beta <= *delta; ++beta){
                        if(beta_vector->at(beta) == 0){
                            continue;
                        }
                        if(i <= beta_vector->at(beta) && j <= beta){
                            flag = true;
                            break;
                        }
                    }
                }
                if(!flag){
                    new_right_index_map->at(r)->insert(j,i);
                }
            }
            r_index->clear();
        }

        *delta = find_core(B);


        for(const auto &l:*isolated_left_vertex_set){
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }

        for(const auto &r:*isolated_right_vertex_set){
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }

        swap(*left_index_map, *new_left_index_map);
        swap(*right_index_map, *new_right_index_map);
    }

    void edge_bipartite_core_maintenance::batch_remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                       const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                       const shared_ptr<uint32_t> &delta,
                                                       const shared_ptr<thread_pool> &pool)
    {
        auto isolated_left_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto isolated_right_vertex_set = make_shared<unordered_set<uint32_t>>();
        B->remove_edge_collection(edge_set, isolated_left_vertex_set, isolated_right_vertex_set);

        auto alpha_vector = make_shared<vector<uint32_t>>(*delta + 1,0);
        for(uint32_t alpha = 1;alpha <= *delta; ++alpha){
            pool->submit_task([=]{
                uint32_t tau_alpha = 0;
                for(const auto &e:*edge_set){
                    auto l = e->get_left_vertex_id();
                    auto r = e->get_right_vertex_id();

                    auto max_value = std::min(left_index_map->at(l)->get_j(alpha), right_index_map->at(r)->get_maximal_j(alpha));
                    if(max_value > tau_alpha){
                        tau_alpha = max_value;
                    }
                }
                alpha_vector->at(alpha) = tau_alpha;
           });
        }

        auto beta_vector = make_shared<vector<uint32_t>>(*delta + 1, 0);
        for(uint32_t beta = 1; beta <= *delta; ++beta){
            pool->submit_task([=]{
                uint32_t tau_beta = 0;
                for(const auto &e:*edge_set){
                    auto l = e->get_left_vertex_id();
                    auto r = e->get_right_vertex_id();

                    auto max_value = std::min(left_index_map->at(l)->get_maximal_i(beta), right_index_map->at(r)->get_i(beta));
                    if(max_value > tau_beta){
                        tau_beta = max_value;
                    }
                }
                beta_vector->at(beta) = tau_beta;
            });
        }
        pool->barrier();

        for(uint32_t alpha = 1;alpha <= *delta; ++alpha){
            if(alpha_vector->at(alpha) == 0){
                continue;
            }
            pool->submit_task([=]{
                batch_update_alpha_tau_alpha_core(B, left_mutex_map, right_mutex_map, new_left_index_map,
                                                  new_right_index_map, alpha, alpha_vector->at(alpha));
            });
        }

        for(uint32_t beta = 1; beta <= *delta; ++beta){
            if(beta_vector->at(beta) == 0){
                continue;
            }
            pool->submit_task([=]{
                batch_update_tau_beta_beta_core(B, left_mutex_map, right_mutex_map, new_left_index_map,
                                                new_right_index_map, beta_vector->at(beta), beta);
            });
        }
        pool->barrier();

        for(const auto &[l, l_index]:*left_index_map){
            for (const auto &[i, j]: *l_index->get_index_map()) {
                bool flag = false;
                for (uint32_t alpha = 1; alpha <= *delta; ++alpha) {
                    if (alpha_vector->at(alpha) == 0) {
                        continue;
                    }
                    if (i <= alpha && j <= alpha_vector->at(alpha)) {
                        flag = true;
                        break;
                    }
                }
                if (!flag) {
                    for (uint32_t beta = 1; beta <= *delta; ++beta) {
                        if (beta_vector->at(beta) == 0) {
                            continue;
                        }
                        if (i <= beta_vector->at(beta) && j <= beta) {
                            flag = true;
                            break;
                        }
                    }
                }
                if (!flag) {
                    new_left_index_map->at(l)->insert(i, j);
                }
            }
            l_index->clear();
        }

        for(const auto &[r, r_index]:*right_index_map){
            for (const auto &[j, i]: *r_index->get_index_map()) {
                bool flag = false;
                for (uint32_t alpha = 1; alpha <= *delta; ++alpha) {
                    if (alpha_vector->at(alpha) == 0) {
                        continue;
                    }
                    if (i <= alpha && j <= alpha_vector->at(alpha)) {
                        flag = true;
                        break;
                    }
                }
                if (!flag) {
                    for (uint32_t beta = 1; beta <= *delta; ++beta) {
                        if (beta_vector->at(beta) == 0) {
                            continue;
                        }
                        if (i <= beta_vector->at(beta) && j <= beta) {
                            flag = true;
                            break;
                        }
                    }
                }
                if (!flag) {
                    new_right_index_map->at(r)->insert(j, i);
                }
            }
            r_index->clear();
        }

        *delta = find_core(B);

        for(const auto &l:*isolated_left_vertex_set){
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }

        for(const auto &r:*isolated_right_vertex_set){
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }

        swap(*left_index_map, *new_left_index_map);
        swap(*right_index_map, *new_right_index_map);
    }

    void edge_bipartite_core_maintenance::batch_update_alpha_phi_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                            uint32_t alpha,
                                                                            uint32_t phi_alpha) {
        auto left_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for(const auto &[l,l_vertex]:*B->get_left_vertex_map()){
            if(l_vertex->get_degree() < alpha){
                evicted_l_set->insert(l);
            }else
            {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            if(r_vertex->get_degree() < phi_alpha){
                evicted_r_set->insert(r);
            } else
            {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);
                left_degree_map->erase(l);

                for(const auto &[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(right_degree_map->count(r) && right_degree_map->at(r) >= phi_alpha){
                        --right_degree_map->at(r);

                        if(right_degree_map->at(r) < phi_alpha){
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= alpha){
                        --left_degree_map->at(l);

                        if(left_degree_map->at(l) < alpha){
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }


        while(!right_degree_map->empty()){
            auto phi = right_degree_map->begin()->second;
            for(const auto &[r, r_degree]:*right_degree_map){
                if(r_degree < phi){
                    phi = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                }else if(r_degree == phi){
                    evicted_r_set->insert(r);
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                for(uint32_t k = 1; k <= phi;++k){
                    new_right_index_map->at(r)->insert(k,alpha);
                }

                for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= alpha){
                        --left_degree_map->at(l);
                        if(left_degree_map->at(l) < alpha){
                            left_degree_map->erase(l);

                            new_left_index_map->at(l)->insert(alpha, phi);

                            for(const auto &[r2,e2]:*B->get_left_vertex(l)->get_edge_map()){
                                if(right_degree_map->count(r2) && right_degree_map->at(r2) > phi){
                                    --right_degree_map->at(r2);
                                    if(right_degree_map->at(r2) <= phi){
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


    void edge_bipartite_core_maintenance::batch_update_alpha_phi_alpha_core(
            const shared_ptr<abstract_bipartite_graph> &B,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
            uint32_t alpha,
            uint32_t phi_alpha) {
        auto left_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for(const auto &[l,l_vertex]:*B->get_left_vertex_map()){
            if(l_vertex->get_degree() < alpha){
                evicted_l_set->insert(l);
            }else
            {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            if(r_vertex->get_degree() < phi_alpha){
                evicted_r_set->insert(r);
            } else
            {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);
                left_degree_map->erase(l);

                for(const auto &[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(right_degree_map->count(r) && right_degree_map->at(r) >= phi_alpha){
                        --right_degree_map->at(r);

                        if(right_degree_map->at(r) < phi_alpha){
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= alpha){
                        --left_degree_map->at(l);

                        if(left_degree_map->at(l) < alpha){
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }


        while(!right_degree_map->empty()){
            auto phi = right_degree_map->begin()->second;
            for(const auto &[r, r_degree]:*right_degree_map){
                if(r_degree < phi){
                    phi = r_degree;

                    evicted_r_set->clear();
                    evicted_r_set->insert(r);
                }else if(r_degree == phi){
                    evicted_r_set->insert(r);
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                right_mutex_map->at(r)->lock();
                for(uint32_t k = 1; k <= phi;++k){
                    new_right_index_map->at(r)->insert(k,alpha);
                }
                right_mutex_map->at(r)->unlock();

                for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= alpha){
                        --left_degree_map->at(l);
                        if(left_degree_map->at(l) < alpha){
                            left_degree_map->erase(l);

                            left_mutex_map->at(l)->lock();
                            new_left_index_map->at(l)->insert(alpha, phi);
                            left_mutex_map->at(l)->unlock();

                            for(const auto &[r2,e2]:*B->get_left_vertex(l)->get_edge_map()){
                                if(right_degree_map->count(r2) && right_degree_map->at(r2) > phi){
                                    --right_degree_map->at(r2);
                                    if(right_degree_map->at(r2) <= phi){
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

    void edge_bipartite_core_maintenance::batch_update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                          uint32_t phi_beta,
                                                                          uint32_t beta) {
        auto left_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for(const auto &[l,l_vertex]:*B->get_left_vertex_map()){
            if(l_vertex->get_degree() < phi_beta){
                evicted_l_set->insert(l);
            } else
            {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            if(r_vertex->get_degree() < beta){
                evicted_r_set->insert(r);
            }else
            {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);
                left_degree_map->erase(l);

                for(const auto &[r, e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(right_degree_map->count(r) && right_degree_map->at(r) >= beta){
                        --right_degree_map->at(r);
                        if(right_degree_map->at(r) < beta){
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                for(const auto &[l, e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= phi_beta){
                        --left_degree_map->at(l);
                        if(left_degree_map->at(l) < phi_beta){
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }


        while(!left_degree_map->empty()){
            uint32_t phi = left_degree_map->begin()->second;
            for(const auto &[l, l_degree]:*left_degree_map){
                if(l_degree < phi){
                    phi = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                }else if(l_degree == phi){
                    evicted_l_set->insert(l);
                }
            }

            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_degree_map->erase(l);
                for (uint32_t k = 1; k <= phi; ++k) {
                    new_left_index_map->at(l)->insert(k, beta);
                }

                for(const auto &[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(right_degree_map->count(r) && right_degree_map->at(r) >= beta){
                        --right_degree_map->at(r);
                        if(right_degree_map->at(r) < beta){
                            right_degree_map->erase(r);

                            new_right_index_map->at(r)->insert(beta,phi);

                            for(const auto &[l2,e2]:*B->get_right_vertex(r)->get_edge_map()){
                                if(left_degree_map->count(l2) && left_degree_map->at(l2) > phi){
                                    --left_degree_map->at(l2);
                                    if(left_degree_map->at(l2) <= phi){
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

    void edge_bipartite_core_maintenance::batch_update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                          uint32_t phi_beta,
                                                                          uint32_t beta) {
        auto left_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for(const auto &[l,l_vertex]:*B->get_left_vertex_map()){
            if(l_vertex->get_degree() < phi_beta){
                evicted_l_set->insert(l);
            } else
            {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            if(r_vertex->get_degree() < beta){
                evicted_r_set->insert(r);
            }else
            {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        }

        while(!evicted_l_set->empty() || !evicted_r_set->empty()){
            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);
                left_degree_map->erase(l);

                for(const auto &[r, e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(right_degree_map->count(r) && right_degree_map->at(r) >= beta){
                        --right_degree_map->at(r);
                        if(right_degree_map->at(r) < beta){
                            evicted_r_set->insert(r);
                        }
                    }
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                for(const auto &[l, e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= phi_beta){
                        --left_degree_map->at(l);
                        if(left_degree_map->at(l) < phi_beta){
                            evicted_l_set->insert(l);
                        }
                    }
                }
            }
        }


        while(!left_degree_map->empty()){
            uint32_t phi = left_degree_map->begin()->second;
            for(const auto &[l, l_degree]:*left_degree_map){
                if(l_degree < phi){
                    phi = l_degree;

                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                }else if(l_degree == phi){
                    evicted_l_set->insert(l);
                }
            }

            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_degree_map->erase(l);

                left_mutex_map->at(l)->lock();
                for (uint32_t k = 1; k <= phi; ++k) {
                    new_left_index_map->at(l)->insert(k, beta);
                }
                left_mutex_map->at(l)->unlock();

                for(const auto &[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(right_degree_map->count(r) && right_degree_map->at(r) >= beta){
                        --right_degree_map->at(r);
                        if(right_degree_map->at(r) < beta){
                            right_degree_map->erase(r);

                            right_mutex_map->at(r)->lock();
                            new_right_index_map->at(r)->insert(beta, phi);
                            right_mutex_map->at(r)->unlock();

                            for(const auto &[l2,e2]:*B->get_right_vertex(r)->get_edge_map()){
                                if(left_degree_map->count(l2) && left_degree_map->at(l2) > phi){
                                    --left_degree_map->at(l2);
                                    if(left_degree_map->at(l2) <= phi){
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

    void edge_bipartite_core_maintenance::batch_update_alpha_tau_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                            uint32_t alpha,
                                                                            uint32_t tau_alpha) {
        auto left_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for(const auto &[l,l_vertex]:*B->get_left_vertex_map()){
            if(l_vertex->get_degree() < alpha){
                evicted_l_set->insert(l);
            }else
            {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            right_degree_map->insert({r, r_vertex->get_degree()});
        }

        while(!evicted_l_set->empty() ){
            auto l = *evicted_l_set->begin();
            evicted_l_set->erase(l);
            left_degree_map->erase(l);

            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (right_degree_map->count(r) && right_degree_map->at(r) >= 1) {
                    --right_degree_map->at(r);

                    if (right_degree_map->at(r) < 1) {
                        right_degree_map->erase(r);

                        for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (left_degree_map->count(l2) && left_degree_map->at(l2) >= alpha) {
                                --left_degree_map->at(l2);

                                if (left_degree_map->at(l2) < alpha) {
                                    evicted_l_set->insert(l2);
                                }
                            }
                        }
                    }
                }
            }
        }

        for(uint32_t tau = 1; tau <=tau_alpha;++tau){
            for(const auto &[r, r_degree]:*right_degree_map){
                if(r_degree == tau){
                    evicted_r_set->insert(r);
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                for(uint32_t k = 1; k <= tau; ++k){
                    new_right_index_map->at(r)->insert(k, alpha);
                }

                for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= alpha){
                        --left_degree_map->at(l);
                        if(left_degree_map->at(l) < alpha){
                            left_degree_map->erase(l);

                            new_left_index_map->at(l)->insert(alpha, tau);

                            for(const auto &[r2,e2]:*B->get_left_vertex(l)->get_edge_map()){
                                if(right_degree_map->count(r2) && right_degree_map->at(r2) > tau){
                                    --right_degree_map->at(r2);
                                    if(right_degree_map->at(r2) <= tau){
                                        evicted_r_set->insert(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for(const auto &[l, l_degree]:*left_degree_map){
            new_left_index_map->at(l)->insert(alpha, tau_alpha);
        }
        for(const auto &[r, r_degree]:*right_degree_map){
            for(uint32_t k = 1; k <= tau_alpha; ++k){
                new_right_index_map->at(r)->insert(k, alpha);
            }
        }
    }

    void edge_bipartite_core_maintenance::batch_update_alpha_tau_alpha_core(
            const shared_ptr<abstract_bipartite_graph> &B,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
            uint32_t alpha,
            uint32_t tau_alpha) {
        auto left_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for(const auto &[l,l_vertex]:*B->get_left_vertex_map()){
            if(l_vertex->get_degree() < alpha){
                evicted_l_set->insert(l);
            }else
            {
                left_degree_map->insert({l, l_vertex->get_degree()});
            }
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            right_degree_map->insert({r, r_vertex->get_degree()});
        }

        while(!evicted_l_set->empty() ){
            auto l = *evicted_l_set->begin();
            evicted_l_set->erase(l);
            left_degree_map->erase(l);

            for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                if (right_degree_map->count(r) && right_degree_map->at(r) >= 1) {
                    --right_degree_map->at(r);

                    if (right_degree_map->at(r) < 1) {
                        right_degree_map->erase(r);

                        for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                            if (left_degree_map->count(l2) && left_degree_map->at(l2) >= alpha) {
                                --left_degree_map->at(l2);

                                if (left_degree_map->at(l2) < alpha) {
                                    evicted_l_set->insert(l2);
                                }
                            }
                        }
                    }
                }
            }
        }

        for(uint32_t tau = 1; tau <=tau_alpha;++tau){
            for(const auto &[r, r_degree]:*right_degree_map){
                if(r_degree == tau){
                    evicted_r_set->insert(r);
                }
            }

            while(!evicted_r_set->empty()){
                auto r = *evicted_r_set->begin();
                evicted_r_set->erase(r);
                right_degree_map->erase(r);

                right_mutex_map->at(r)->lock();
                for(uint32_t k = 1; k <= tau; ++k){
                    new_right_index_map->at(r)->insert(k, alpha);
                }
                right_mutex_map->at(r)->unlock();

                for(const auto &[l,e]:*B->get_right_vertex(r)->get_edge_map()){
                    if(left_degree_map->count(l) && left_degree_map->at(l) >= alpha){
                        --left_degree_map->at(l);
                        if(left_degree_map->at(l) < alpha){
                            left_degree_map->erase(l);

                            left_mutex_map->at(l)->lock();
                            new_left_index_map->at(l)->insert(alpha, tau);
                            left_mutex_map->at(l)->unlock();

                            for(const auto &[r2,e2]:*B->get_left_vertex(l)->get_edge_map()){
                                if(right_degree_map->count(r2) && right_degree_map->at(r2) > tau){
                                    --right_degree_map->at(r2);
                                    if(right_degree_map->at(r2) == tau){
                                        evicted_r_set->insert(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for(const auto &[l, l_degree]:*left_degree_map){
            left_mutex_map->at(l)->lock();
            new_left_index_map->at(l)->insert(alpha, tau_alpha);
            left_mutex_map->at(l)->unlock();
        }

        for(const auto &[r, r_degree]:*right_degree_map){
            right_mutex_map->at(r)->lock();
            for(uint32_t k = 1; k <= tau_alpha; ++k){
                new_right_index_map->at(r)->insert(k, alpha);
            }
            right_mutex_map->at(r)->unlock();
        }
    }

    void edge_bipartite_core_maintenance::batch_update_tau_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                          uint32_t tau_beta,
                                                                          uint32_t beta) {
        auto left_degree_map = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for(const auto &[l,l_vertex]:*B->get_left_vertex_map()){
            left_degree_map->insert({l, l_vertex->get_degree()});
        }

        for(const auto &[r, r_vertex]:*B->get_right_vertex_map()){
            if(r_vertex->get_degree() < beta){
                evicted_r_set->insert(r);
            }else
            {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        }

        while(!evicted_r_set->empty()){
            auto r = *evicted_r_set->begin();
            evicted_r_set->erase(r);
            right_degree_map->erase(r);

            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if (left_degree_map->count(l) && left_degree_map->at(l) >= 1) {
                    --left_degree_map->at(l);

                    if (left_degree_map->at(l) < 1) {
                        left_degree_map->erase(l);

                        for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (right_degree_map->count(r2) && right_degree_map->at(r2) >= beta) {
                                --right_degree_map->at(r2);
                                if (right_degree_map->at(r2) < beta) {
                                    evicted_r_set->insert(r2);
                                }
                            }
                        }
                    }
                }
            }
        }


        for(uint32_t tau =1; tau <=tau_beta;++tau){
            for(const auto &[l, l_degree]:*left_degree_map){
                if(l_degree == tau){
                    evicted_l_set->insert(l);
                }
            }

            while(!evicted_l_set->empty()){
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_degree_map->erase(l);

                for(uint32_t k = 1; k <= tau; ++k) {
                    left_index_map->at(l)->insert(k, beta);
                }


                for(const auto &[r,e]:*B->get_left_vertex(l)->get_edge_map()){
                    if(right_degree_map->count(r) && right_degree_map->at(r) >= beta){
                        --right_degree_map->at(r);
                        if(right_degree_map->at(r) < beta){
                            right_degree_map->erase(r);

                            right_index_map->at(r)->insert(beta, tau);

                            for(const auto &[l2,e2]:*B->get_right_vertex(r)->get_edge_map()){
                                if(left_degree_map->count(l2) && left_degree_map->at(l2) > tau){
                                    --left_degree_map->at(l2);
                                    if(left_degree_map->at(l2) == tau){
                                        evicted_l_set->insert(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for(const auto &[l, l_degree]:*left_degree_map){
            for(uint32_t k = 1; k <= tau_beta; ++k) {
                left_index_map->at(l)->insert(k, beta);
            }
        }

        for(const auto &[r, r_degree]:*right_degree_map){
            right_index_map->at(r)->insert(beta, tau_beta);
        }
    }

    void edge_bipartite_core_maintenance::batch_update_tau_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                          uint32_t tau_beta,
                                                                          uint32_t beta) {
        auto left_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        for (const auto &[l, l_vertex]: *B->get_left_vertex_map()) {
            left_degree_map->insert({l, l_vertex->get_degree()});
        }

        for (const auto &[r, r_vertex]: *B->get_right_vertex_map()) {
            if (r_vertex->get_degree() < beta) {
                evicted_r_set->insert(r);
            } else {
                right_degree_map->insert({r, r_vertex->get_degree()});
            }
        }

        while (!evicted_r_set->empty()) {
            auto r = *evicted_r_set->begin();
            evicted_r_set->erase(r);
            right_degree_map->erase(r);

            for (const auto &[l, e]: *B->get_right_vertex(r)->get_edge_map()) {
                if (left_degree_map->count(l) && left_degree_map->at(l) >= 1) {
                    --left_degree_map->at(l);

                    if (left_degree_map->at(l) < 1) {
                        left_degree_map->erase(l);

                        for (const auto &[r2, e2]: *B->get_left_vertex(l)->get_edge_map()) {
                            if (right_degree_map->count(r2) && right_degree_map->at(r2) >= beta) {
                                --right_degree_map->at(r2);
                                if (right_degree_map->at(r2) < beta) {
                                    evicted_r_set->insert(r2);
                                }
                            }
                        }
                    }
                }
            }
        }


        for (uint32_t tau = 1; tau <= tau_beta; ++tau) {
            for (const auto &[l, l_degree]: *left_degree_map) {
                if (l_degree == tau) {
                    evicted_l_set->insert(l);
                }
            }

            while (!evicted_l_set->empty()) {
                auto l = *evicted_l_set->begin();
                evicted_l_set->erase(l);

                left_degree_map->erase(l);

                left_mutex_map->at(l)->lock();
                for (uint32_t k = 1; k <= tau; ++k) {
                    new_left_index_map->at(l)->insert(k, beta);
                }
                left_mutex_map->at(l)->unlock();


                for (const auto &[r, e]: *B->get_left_vertex(l)->get_edge_map()) {
                    if (right_degree_map->count(r) && right_degree_map->at(r) >= beta) {
                        --right_degree_map->at(r);
                        if (right_degree_map->at(r) < beta) {
                            right_degree_map->erase(r);

                            right_mutex_map->at(r)->lock();
                            new_right_index_map->at(r)->insert(beta, tau);
                            right_mutex_map->at(r)->unlock();

                            for (const auto &[l2, e2]: *B->get_right_vertex(r)->get_edge_map()) {
                                if (left_degree_map->count(l2) && left_degree_map->at(l2) > tau) {
                                    --left_degree_map->at(l2);
                                    if (left_degree_map->at(l2) == tau) {
                                        evicted_l_set->insert(l2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (const auto &[l, l_degree]: *left_degree_map) {
            left_mutex_map->at(l)->lock();
            for (uint32_t k = 1; k <= tau_beta; ++k) {
                new_left_index_map->at(l)->insert(k, beta);
            }
            left_mutex_map->at(l)->unlock();
        }

        for (const auto &[r, r_degree]: *right_degree_map) {
            right_mutex_map->at(r)->lock();
            new_right_index_map->at(r)->insert(beta, tau_beta);
            right_mutex_map->at(r)->unlock();
        }
    }

    uint32_t edge_bipartite_core_maintenance::find_core(const shared_ptr<abstract_bipartite_graph> &B) {
        auto left_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto right_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for(const auto&[l,l_vertex]:*B->get_left_vertex_map()){
            left_degree_map->insert({l, l_vertex->get_degree()});
        }
        for(const auto&[r,r_vertex]:*B->get_right_vertex_map()){
            right_degree_map->insert({r, r_vertex->get_degree()});
        }

        auto evicted_l_set = make_shared<unordered_set<uint32_t>>();
        auto evicted_r_set = make_shared<unordered_set<uint32_t>>();

        uint32_t max_delta = 1;
        while(!left_degree_map->empty() && !right_degree_map->empty())
        {
            uint32_t  delta = left_degree_map->begin()->second;
            for (const auto& [l, l_degree]:*left_degree_map) {
                if (l_degree < delta) {
                    delta = l_degree;
                    evicted_l_set->clear();
                    evicted_l_set->insert(l);
                } else if(l_degree == delta){
                    evicted_l_set->insert(l);
                }
            }

            for (const auto& [r,r_degree]:*right_degree_map) {
                if (r_degree < delta) {
                    delta = r_degree;

                    evicted_l_set->clear();
                    evicted_r_set->clear();

                    evicted_r_set->insert(r);
                }else if(r_degree == delta){
                    evicted_r_set->insert(r);
                }
            }

            max_delta = delta;

            while (!evicted_l_set->empty() || !evicted_r_set->empty()) {
                while (!evicted_l_set->empty()) {
                    auto l = *evicted_l_set->begin();
                    evicted_l_set->erase(l);
                    left_degree_map->erase(l);

                    for (const auto&[r,e]:*B->get_left_vertex(l)->get_edge_map()) {
                        if(right_degree_map->count(r) && right_degree_map->at(r) > delta){
                            --right_degree_map->at(r);
                            if (right_degree_map->at(r) == delta) {
                                evicted_r_set->insert(r);
                            }
                        }
                    }
                }

                while (!evicted_r_set->empty()) {
                    auto r = *evicted_r_set->begin();
                    evicted_r_set->erase(r);
                    right_degree_map->erase(r);

                    for (const auto& [l,e]:*B->get_right_vertex(r)->get_edge_map()) {
                        if(left_degree_map->count(l) && left_degree_map->at(l) > delta){
                            --left_degree_map->at(l);
                            if (left_degree_map->at(l) == delta) {
                                evicted_l_set->insert(l);
                            }
                        }
                    }
                }
            }
        }
        return max_delta;
    }

    uint32_t edge_bipartite_core_maintenance::get_maximal_b_alpha(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                  uint32_t l,
                                                                  uint32_t alpha)
    {
        auto count_map = make_shared<map<uint32_t,uint32_t>>();
        auto l_vertex = B->get_left_vertex(l);
        if(!l_vertex){
            return 0;
        }
        for(const auto&[r,e]:*l_vertex->get_edge_map()){
            auto r_index_map = right_index_map->at(r)->get_index_map();
            for (const auto&[j, i]:*r_index_map) {
                if (i >= alpha) {
                    if (!count_map->count(j)) {
                        count_map->insert({j, 0});
                    }
                    ++count_map->at(j);
                }
            }
        }

        uint32_t b_alpha = 0;
        for(auto iter = count_map->rbegin();iter!=count_map->rend();++iter){
            auto x = iter->first;
            auto count = iter->second;
            if(count >= alpha){
                b_alpha = x;
                break;
            }
        }

        return b_alpha;
    }

    uint32_t edge_bipartite_core_maintenance::get_maximal_b_beta(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                 uint32_t r,
                                                                 uint32_t beta)
    {
        auto count_map = make_shared<map<uint32_t, uint32_t>>();
        auto r_vertex = B->get_right_vertex(r);
        if(!r_vertex){
            return 0;
        }
        for (const auto&[l, e]:*r_vertex->get_edge_map()) {
            auto l_map = left_index_map->at(l)->get_index_map();
            for (const auto&[i, j]:*l_map) {
                if (j >= beta) {
                    if (!count_map->count(i)) {
                        count_map->insert({i, 0});
                    }
                    ++count_map->at(i);
                }
            }
        }

        uint32_t b_beta = 0;
        for(auto iter = count_map->rbegin();iter!=count_map->rend();++iter){
            auto x = iter->first;
            auto count = iter->second;
            if(count >= beta){
                b_beta = x;
                break;
            }
        }
        return b_beta;
    }

    void edge_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map) {
        for(const auto &[l, l_vertex]:*B->get_left_vertex_map()){
            new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
        }

        for(const auto &[r, r_index]:*B->get_right_vertex_map()) {
            new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
        }
    }
    
    void edge_bipartite_core_maintenance::init(const shared_ptr<abstract_bipartite_graph> &B,
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
        for(const auto &[r, r_index]:*B->get_right_vertex_map()){
            pool->submit_task([=]{
                new_right_index_map->at(r) = make_shared<bipartite_core_right_store_index>();
                right_mutex_map->at(r) = make_shared<mutex>();
            });
        }
        pool->barrier();
    }

    void edge_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<abstract_bipartite_edge> &e,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                 const shared_ptr<uint32_t>& delta) {
        {
            B->insert_edge(e);

            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(!left_index_map->count(l))
            {
                left_index_map->insert({l,make_shared<bipartite_core_left_store_index>()});
                left_index_map->at(l)->insert(1, 1);

                new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
            }
            if(!right_index_map->count(r)){
                right_index_map->insert({r,make_shared<bipartite_core_right_store_index>()});
                right_index_map->at(r)->insert(1, 1);

                new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
            }
        }

        *delta = find_core(B);

        for (uint32_t alpha = 1; alpha <= *delta; ++alpha) {
            update_alpha_phi_alpha_core(B, left_index_map, right_index_map,
                                        new_left_index_map, new_right_index_map, e, alpha);
        }

        for(uint32_t beta = 1; beta <= *delta; ++beta){
            update_phi_beta_beta_core(B, left_index_map, right_index_map,
                                      new_left_index_map,
                                      new_right_index_map, e, beta);
        }

        for (const auto &[l1, l1_index]: *new_left_index_map) {
            for (const auto &[i, j]: *l1_index->get_index_map()) {
                left_index_map->at(l1)->insert(i, j);
            }
            new_right_index_map->at(l1)->clear();
        }

        for (const auto &[r1, r1_index]: *new_right_index_map) {
            for (const auto &[j, i]: *r1_index->get_index_map()) {
                right_index_map->at(r1)->insert(j, i);
            }
            new_right_index_map->at(r1)->clear();
        }
    }

    void edge_bipartite_core_maintenance::insert(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                 const shared_ptr<abstract_bipartite_edge> &e,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                 const shared_ptr<uint32_t> &delta,
                                                 const shared_ptr<thread_pool> &pool) {
        {
            B->insert_edge(e);

            auto l = e->get_left_vertex_id();
            auto r = e->get_right_vertex_id();

            if(!left_index_map->count(l))
            {
                left_index_map->insert({l,make_shared<bipartite_core_left_store_index>()});
                left_index_map->at(l)->insert(1, 1);

                new_left_index_map->insert({l, make_shared<bipartite_core_left_store_index>()});
                left_mutex_map->insert({l, make_shared<mutex>()});
            }
            if(!right_index_map->count(r)){
                right_index_map->insert({r,make_shared<bipartite_core_right_store_index>()});
                right_index_map->at(r)->insert(1, 1);

                new_right_index_map->insert({r, make_shared<bipartite_core_right_store_index>()});
                right_mutex_map->insert({r, make_shared<mutex>()});
            }
        }
        *delta = find_core(B);

        for(auto alpha = 1;alpha <= *delta;++alpha){
            pool->submit_task([=]{
                update_alpha_phi_alpha_core(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map,
                                            new_left_index_map, new_right_index_map, e, alpha);
            });
        }

        for(auto beta = 1;beta <= *delta;++beta){
            pool->submit_task([=]{
                update_phi_beta_beta_core(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map,
                                          new_left_index_map,
                                          new_right_index_map, e, beta);
            });
        }
        pool->barrier();

        for (const auto &[l, l_index]: *new_left_index_map) {
            for (const auto &[i, j]: *l_index->get_index_map()) {
                left_index_map->at(l)->insert(i, j);
            }
            new_left_index_map->at(l)->clear();
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            for (const auto &[j, i]: *r_index->get_index_map()) {
                right_index_map->at(r)->insert(j, i);
            }
            new_right_index_map->at(r)->clear();
        }
    }

    void edge_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<abstract_bipartite_edge> &e,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                 const shared_ptr<uint32_t>& delta)
    {
        B->remove_edge(e);

        for (auto alpha = 1; alpha <= *delta; ++alpha) {
            update_alpha_tau_alpha_core(B, left_index_map, right_index_map, new_left_index_map, new_right_index_map, e,
                                        alpha);
        }

        for (auto beta = 1; beta <= *delta; ++beta) {
            update_beta_tau_beta_core(B, left_index_map, right_index_map, new_left_index_map, new_right_index_map, e,
                                      beta);
        }

        for (const auto &[l, l_index]: *new_left_index_map) {
            for (const auto &[i, j]: *l_index->get_index_map()) {
                left_index_map->at(l)->remove(i, j);
                if(left_index_map->at(l)->get_j(i) == 0){
                    left_index_map->at(l)->remove(i);
                }
            }
            new_left_index_map->at(l)->clear();
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            for (const auto &[j, i]: *r_index->get_index_map()) {
                right_index_map->at(r)->remove(j, i);
                if(right_index_map->at(r)->get_i(j) == 0){
                    right_index_map->at(r)->remove(j);
                }
            }
            new_right_index_map->at(r)->clear();
        }


        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        if(!B->get_left_vertex(l)) {
            left_index_map->erase(l);
            new_left_index_map->erase(l);
        }
        if(!B->get_right_vertex(r)){
            right_index_map->erase(r);
            new_right_index_map->erase(r);
        }

        *delta = find_core(B);
    }

    void edge_bipartite_core_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &B,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                 const shared_ptr<abstract_bipartite_edge> &e,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                 const shared_ptr<uint32_t> &delta,
                                                 const shared_ptr<thread_pool> &pool)
    {
        B->remove_edge(e);

        for (auto alpha = 1; alpha <= *delta; ++alpha) {
            pool->submit_task([=]{
                update_alpha_tau_alpha_core(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map,
                                            new_left_index_map, new_right_index_map, e,
                                            alpha);
            });
        }

        for (auto beta = 1; beta <= *delta; ++beta) {
            pool->submit_task([=]{
                update_beta_tau_beta_core(B, left_mutex_map, right_mutex_map, left_index_map, right_index_map,
                                          new_left_index_map, new_right_index_map, e,
                                          beta);
            });
        }
        pool->barrier();


        for (const auto &[l, l_index]: *new_left_index_map) {
            for (const auto &[i, j]: *l_index->get_index_map()) {
                left_index_map->at(l)->remove(i, j);
                if(left_index_map->at(l)->get_j(i) == 0){
                    left_index_map->at(l)->remove(i);
                }
            }
            new_left_index_map->at(l)->clear();
        }

        for (const auto &[r, r_index]: *new_right_index_map) {
            for (const auto &[j, i]: *r_index->get_index_map()) {
                right_index_map->at(r)->remove(j, i);
                if(right_index_map->at(r)->get_i(j) == 0){
                    right_index_map->at(r)->remove(j);
                }
            }
            new_right_index_map->at(r)->clear();
        }

        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        if(!B->get_left_vertex(l)) {
            left_index_map->erase(l);
            new_left_index_map->erase(l);
            left_mutex_map->erase(l);
        }
        if(!B->get_right_vertex(r)){
            right_index_map->erase(r);
            new_right_index_map->erase(r);
            right_mutex_map->erase(r);
        }

        *delta = find_core(B);
    }



    void edge_bipartite_core_maintenance::remove_left_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cl,
                                                                 const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cr,
                                                                 uint32_t l,
                                                                 uint32_t alpha,
                                                                 uint32_t beta)
    {
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        visited_set->insert(l);

        while(!visited_set->empty()){
            auto l1 = *visited_set->begin();
            visited_set->erase(l1);
            Cl->erase(l1);

            for(const auto&[r1,e1]:* B->get_left_vertex(l1)->get_edge_map()){
                if(Cr->count(r1) && Cr->at(r1) > 0){
                    --Cr->at(r1);
                    if(Cr->at(r1) < beta){
                        Cr->erase(r1);

                        for(const auto&[l2,e2]:*B->get_right_vertex(r1)->get_edge_map()){
                            if(Cl->count(l2) && Cl->at(l2) > 0) {
                                --Cl->at(l2);
                                if(Cl->at(l2) < alpha){
                                    visited_set->insert(l2);
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    void edge_bipartite_core_maintenance::remove_right_candidates(const shared_ptr<abstract_bipartite_graph> &B,
                                                                  const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cl,
                                                                  const shared_ptr<unordered_map<uint32_t,uint32_t>> &Cr,
                                                                  uint32_t r,
                                                                  uint32_t alpha,
                                                                  uint32_t beta)
    {
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        visited_set->insert(r);

        while(!visited_set->empty()){
            auto r1 = *visited_set->begin();
            visited_set->erase(r1);
            Cr->erase(r1);


            for(const auto&[l1,e]:*B->get_right_vertex(r1)->get_edge_map()){
                if(Cl->count(l1) && Cl->at(l1) > 0) {
                    --Cl->at(l1);
                    if(Cl->at(l1) < alpha){
                        Cl->erase(l1);

                        for(const auto &[r2,e2]:*B->get_left_vertex(l1)->get_edge_map()){
                            if(Cr->count(r2) && Cr->at(r2) > 0){
                                --Cr->at(r2);

                                if(Cr->at(r2) < beta){
                                    visited_set->insert(r2);
                                }
                            }
                        }
                    }
                }
            }
        }


    }

    void edge_bipartite_core_maintenance::update_alpha_phi_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                      const shared_ptr<abstract_bipartite_edge>& e,
                                                                      uint32_t alpha){

        auto right_vertex_beta_max_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_alpha = get_maximal_b_alpha(B, right_index_map, l, alpha);

        if(!right_vertex_beta_max_map->count(r)){
            right_vertex_beta_max_map->insert({r, right_index_map->at(r)->get_maximal_j(alpha)});
        }
        auto phi_alpha = min(b_alpha, right_vertex_beta_max_map->at(r));
        if(phi_alpha == 0){
            return;
        }

        new_left_index_map->at(l)->insert(alpha, b_alpha);

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto Cl = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto Cr = make_shared<unordered_map<uint32_t,uint32_t>>();

        auto l_value = b_alpha > left_index_map->at(l)->get_j(alpha) ? b_alpha:left_index_map->at(l)->get_j(alpha);

        if(l_value <= right_vertex_beta_max_map->at(r))
        {
            Sl->insert(l);
        }else
        {
            Sr->insert(r);
        }
        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                Tl->insert(l1);
                Cl->insert({l1, 0});

                auto l1_vertex = B->get_left_vertex(l1);
                for(const auto&[r1,e1]:*l1_vertex->get_edge_map())
                {
                    if(!right_vertex_beta_max_map->count(r1)){
                        right_vertex_beta_max_map->insert({r1, right_index_map->at(r1)->get_maximal_j(alpha)});
                    }
                    auto beta_max = right_vertex_beta_max_map->at(r1);
                    if(beta_max == 0){
                        continue;
                    }
                    if(beta_max > phi_alpha || Cr->count(r1))
                    {
                        ++Cl->at(l1);
                    }else if(!Tr->count(r1) && beta_max == phi_alpha){
                        ++Cl->at(l1);
                        Sr->insert(r1);
                    }
                }
                if(Cl->at(l1) < alpha)
                {
                    remove_left_candidates(B, Cl, Cr, l1, alpha, phi_alpha + 1);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                Tr->insert(r1);
                Cr->insert({r1, 0});

                for(const auto& [l1,e1]:*B->get_right_vertex(r1)->get_edge_map()){

                    auto beta_max = l == l1 ? l_value:left_index_map->at(l1)->get_j(alpha);
                    if(beta_max == 0){
                        continue;
                    }
                    if(beta_max > phi_alpha || Cl->count(l1)){
                        ++Cr->at(r1);
                    }else if(!Tl->count(l1) && beta_max == phi_alpha){
                        ++Cr->at(r1);
                        Sl->insert(l1);
                    }
                }

                if(Cr->at(r1) < phi_alpha + 1){
                    remove_right_candidates(B, Cl, Cr, r1, alpha, phi_alpha + 1);
                }
            }
        }

        for(const auto&[r1,r1_degree]:*Cr)
        {
            for(uint32_t j = 1;j <= phi_alpha+1;++j) {
                new_right_index_map->at(r1)->insert(j, alpha);
            }
        }

        for(const auto&[l1,degree]:*Cl){
            new_left_index_map->at(l1)->insert(alpha, phi_alpha + 1);
        }
    }

    void edge_bipartite_core_maintenance::update_alpha_phi_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                      const shared_ptr<abstract_bipartite_edge> &e,
                                                                      uint32_t alpha){
        auto right_vertex_beta_max_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_alpha = get_maximal_b_alpha(B, right_index_map, l, alpha);

        if(!right_vertex_beta_max_map->count(r)){
            right_vertex_beta_max_map->insert({r, right_index_map->at(r)->get_maximal_j(alpha)});
        }
        auto phi_alpha = min(b_alpha, right_vertex_beta_max_map->at(r));
        if(phi_alpha == 0){
            return;
        }

        left_mutex_map->at(l)->lock();
        new_left_index_map->at(l)->insert(alpha, b_alpha);
        left_mutex_map->at(l)->unlock();

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto Cl = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto Cr = make_shared<unordered_map<uint32_t,uint32_t>>();

        auto l_value = b_alpha > left_index_map->at(l)->get_j(alpha) ? b_alpha:left_index_map->at(l)->get_j(alpha);

        if(l_value <= right_vertex_beta_max_map->at(r))
        {
            Sl->insert(l);
        }else
        {
            Sr->insert(r);
        }

        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                Tl->insert(l1);
                Cl->insert({l1, 0});

                auto l1_vertex = B->get_left_vertex(l1);
                for(const auto&[r1,e1]:*l1_vertex->get_edge_map())
                {
                    if(!right_vertex_beta_max_map->count(r1)){
                        right_vertex_beta_max_map->insert({r1, right_index_map->at(r1)->get_maximal_j(alpha)});
                    }
                    auto beta_max = right_vertex_beta_max_map->at(r1);
                    if(beta_max == 0){
                        continue;
                    }
                    if(beta_max > phi_alpha || Cr->count(r1))
                    {
                        ++Cl->at(l1);
                    }else if(!Tr->count(r1) && beta_max == phi_alpha){
                        ++Cl->at(l1);
                        Sr->insert(r1);
                    }
                }
                if(Cl->at(l1) < alpha)
                {
                    remove_left_candidates(B, Cl, Cr ,l1, alpha, phi_alpha + 1);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                Tr->insert(r1);
                Cr->insert({r1, 0});

                for(const auto& [l1,e1]:*B->get_right_vertex(r1)->get_edge_map()){

                    auto beta_max = l == l1 ? l_value:left_index_map->at(l1)->get_j(alpha);
                    if(beta_max == 0){
                        continue;
                    }
                    if(beta_max > phi_alpha || Cl->count(l1)){
                        ++Cr->at(r1);
                    }else if(!Tl->count(l1) && beta_max == phi_alpha){
                        ++Cr->at(r1);
                        Sl->insert(l1);
                    }
                }

                if(Cr->at(r1) < phi_alpha + 1){
                    remove_right_candidates(B, Cl, Cr,r1, alpha, phi_alpha + 1);
                }
            }
        }


        for(const auto&[r1,r1_degree]:*Cr)
        {
            right_mutex_map->at(r1)->lock();
            for(uint32_t j = 1;j <= phi_alpha+1;++j) {
                new_right_index_map->at(r1)->insert(j, alpha);
            }
            right_mutex_map->at(r1)->unlock();
        }

        for(const auto&[l1,degree]:*Cl){
            left_mutex_map->at(l1)->lock();
            new_left_index_map->at(l1)->insert(alpha, phi_alpha + 1);
            left_mutex_map->at(l1)->unlock();
        }
    }

    void edge_bipartite_core_maintenance::update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph>& B,
                                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t,shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                    const shared_ptr<abstract_bipartite_edge>& e,
                                                                    uint32_t beta){
        auto left_vertex_alpha_max_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_beta = get_maximal_b_beta(B, left_index_map, r, beta);
        if(!left_vertex_alpha_max_map->count(l)){
            left_vertex_alpha_max_map->insert({l,left_index_map->at(l)->get_maximal_i(beta)});
        }
        auto phi_beta = min(b_beta, left_vertex_alpha_max_map->at(l));
        if(phi_beta == 0){
            return;
        }
        /**
         * @brief the pair of the right index map is reverse
         */

        new_right_index_map->at(r)->insert(beta, b_beta);

        uint32_t r_value = b_beta > right_index_map->at(r)->get_i(beta) ? b_beta:right_index_map->at(r)->get_i(beta);
        //right_index_map->at(r)->insert(beta, b_beta);

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto Cl = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto Cr = make_shared<unordered_map<uint32_t,uint32_t>>();

        //right_index_map->at(r)->get_i(beta)
        if(r_value <= left_vertex_alpha_max_map->at(l))
        {
            Sr->insert(r);
        }else
        {
            Sl->insert(l);
        }

        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                Tl->insert(l1);
                Cl->insert({l1, 0});

                for(const auto&[r1,e1]:*B->get_left_vertex(l1)->get_edge_map())
                {

                    //right_index_map->at(r1)->get_i(beta)
                    auto alpha_max = r1 == r ? r_value:right_index_map->at(r1)->get_i(beta);
                    if(alpha_max == 0){
                        continue;
                    }
                    if(alpha_max > phi_beta || Cr->count(r1))
                    {
                        ++Cl->at(l1);
                    }else if(!Tr->count(r1) && alpha_max == phi_beta){
                        ++Cl->at(l1);
                        Sr->insert(r1);
                    }
                }
                if(Cl->at(l1) < phi_beta + 1)
                {
                    remove_left_candidates(B, Cl, Cr, l1, phi_beta + 1, beta);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                Tr->insert(r1);
                Cr->insert({r1, 0});

                for(const auto& [l1,e1]:*B->get_right_vertex(r1)->get_edge_map()){

                    if(!left_vertex_alpha_max_map->count(l1)){
                        left_vertex_alpha_max_map->insert({l1,left_index_map->at(l1)->get_maximal_i(beta)});
                    }
                    auto alpha_max = left_vertex_alpha_max_map->at(l1);
                    if(alpha_max == 0){
                        continue;
                    }
                    if(alpha_max > phi_beta || Cl->count(l1)){
                        ++Cr->at(r1);
                    }else if(!Tl->count(l1) && alpha_max == phi_beta){
                        ++Cr->at(r1);
                        Sl->insert(l1);
                    }
                }

                if(Cr->at(r1) < beta){
                    remove_right_candidates(B, Cl, Cr,  r1, phi_beta + 1, beta);
                }
            }
        }

        for(const auto&[l1,l1_degree]:*Cl)
        {
            for(uint32_t i = 1;i <= phi_beta+1;++i){
                new_left_index_map->at(l1)->insert(i, beta);
            }
        }

        for(const auto&[r1,r1_degree]:*Cr){
            new_right_index_map->at(r1)->insert(beta, phi_beta+1);
        }
    }



    void edge_bipartite_core_maintenance::update_phi_beta_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                    const shared_ptr<abstract_bipartite_edge> &e,
                                                                    uint32_t beta){
        auto left_vertex_alpha_max_map = make_shared<unordered_map<uint32_t,uint32_t>>();

        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_beta = get_maximal_b_beta(B, left_index_map, r, beta);
        if(!left_vertex_alpha_max_map->count(l)){
            left_vertex_alpha_max_map->insert({l,left_index_map->at(l)->get_maximal_i(beta)});
        }
        auto phi_beta = min(b_beta, left_vertex_alpha_max_map->at(l));
        if(phi_beta == 0){
            return;
        }
        /**
         * @brief the pair of the right index map is reverse
         */
        right_mutex_map->at(r)->lock();
        new_right_index_map->at(r)->insert(beta, b_beta);
        right_mutex_map->at(r)->unlock();

        uint32_t r_value = b_beta > right_index_map->at(r)->get_i(beta) ? b_beta:right_index_map->at(r)->get_i(beta);
        //right_index_map->at(r)->insert(beta, b_beta);

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto Cl = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto Cr = make_shared<unordered_map<uint32_t,uint32_t>>();

        //right_index_map->at(r)->get_i(beta)
        if(r_value <= left_vertex_alpha_max_map->at(l))
        {
            Sr->insert(r);
        }else
        {
            Sl->insert(l);
        }

        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                Tl->insert(l1);
                Cl->insert({l1, 0});

                for(const auto&[r1,e1]:*B->get_left_vertex(l1)->get_edge_map())
                {
                    //right_index_map->at(r1)->get_i(beta)
                    auto alpha_max = r1 == r ? r_value:right_index_map->at(r1)->get_i(beta);
                    if(alpha_max == 0){
                        continue;
                    }
                    if(alpha_max > phi_beta || Cr->count(r1))
                    {
                        ++Cl->at(l1);
                    }else if(!Tr->count(r1) && alpha_max == phi_beta){
                        ++Cl->at(l1);
                        Sr->insert(r1);
                    }
                }
                if(Cl->at(l1) < phi_beta + 1)
                {
                    remove_left_candidates(B, Cl, Cr, l1, phi_beta + 1, beta);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                Tr->insert(r1);
                Cr->insert({r1, 0});

                for(const auto& [l1,e1]:*B->get_right_vertex(r1)->get_edge_map()){

                    if(!left_vertex_alpha_max_map->count(l1)){
                        left_vertex_alpha_max_map->insert({l1,left_index_map->at(l1)->get_maximal_i(beta)});
                    }
                    auto alpha_max = left_vertex_alpha_max_map->at(l1);
                    if(alpha_max == 0){
                        continue;
                    }
                    if(alpha_max > phi_beta || Cl->count(l1)){
                        ++Cr->at(r1);
                    }else if(!Tl->count(l1) && alpha_max == phi_beta){
                        ++Cr->at(r1);
                        Sl->insert(l1);
                    }
                }

                if(Cr->at(r1) < beta){
                    remove_right_candidates(B, Cl, Cr, r1, phi_beta + 1, beta);
                }
            }
        }

        for(const auto&[l1,l1_degree]:*Cl)
        {
            left_mutex_map->at(l1)->lock();
            for(uint32_t i = 1;i <= phi_beta+1;++i){
                new_left_index_map->at(l1)->insert(i, beta);
            }
            left_mutex_map->at(l1)->unlock();
        }

        for(const auto&[r1,r1_degree]:*Cr){
            right_mutex_map->at(r1)->lock();
            new_right_index_map->at(r1)->insert(beta, phi_beta + 1);
            right_mutex_map->at(r1)->unlock();
        }
    }

    void edge_bipartite_core_maintenance::update_alpha_tau_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                      const shared_ptr<abstract_bipartite_edge>& e,
                                                                      uint32_t alpha)
    {
        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_alpha = get_maximal_b_alpha(B, right_index_map, l, alpha);
        auto tau_alpha = min(left_index_map->at(l)->get_j(alpha),right_index_map->at(r)->get_maximal_j(alpha));
        if(tau_alpha == 0)
        {
            return;
        }

        auto Cl = make_shared<unordered_set<uint32_t>>();
        auto Cr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto sup_l = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto sup_r = make_shared<unordered_map<uint32_t,uint32_t>>();

        if(left_index_map->at(l)->get_j(alpha) == tau_alpha)
        {
            Sl->insert(l);
        }

        if(right_index_map->at(r)->get_maximal_j(alpha) == tau_alpha)
        {
            Sr->insert(r);
        }

        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                if(!B->get_left_vertex(l1)){
                    continue;
                }

                sup_l->insert({l1, 0});
                Tl->insert(l1);

                for(const auto&[r1,e1]:*B->get_left_vertex(l1)->get_edge_map())
                {
                    auto beta_max = right_index_map->at(r1)->get_maximal_j(alpha);
                    if(beta_max >= tau_alpha && ! Cr->count(r1))
                    {
                        ++sup_l->at(l1);
                    }
                }
                if(sup_l->at(l1) < alpha){
                    add_left_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, l1, alpha, tau_alpha);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                sup_r->insert({r1, 0});
                Tr->insert(r1);

                auto r1_vertex = B->get_right_vertex(r1);
                if(!r1_vertex){
                    continue;
                }
                for(const auto& [l1,e1]:*r1_vertex->get_edge_map()){
                    auto beta_max = left_index_map->at(l1)->get_j(alpha);
                    if(beta_max >= tau_alpha && ! Cl->count(l1)){
                        ++sup_r->at(r1);
                    }
                }
                if(sup_r->at(r1) < tau_alpha){
                    add_right_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, r1, alpha, tau_alpha);
                }
            }
        }

        for(const auto&l1:*Cl) {
            new_left_index_map->at(l1)->remove(alpha, tau_alpha - 1);
        }

        for(const auto& r1:*Cr){
            uint32_t r1_degree = B->get_right_vertex(r1) ? B->get_right_vertex(r1)->get_degree():0;
            if(r1 == r){
                ++r1_degree;
            }
            for(uint32_t j = tau_alpha; j<=r1_degree;++j){
                new_right_index_map->at(r1)->remove(j, alpha - 1);
            }
        }
        //B->remove_edge(e);

        new_left_index_map->at(l)->remove(alpha, b_alpha);
    }

    void edge_bipartite_core_maintenance::update_alpha_tau_alpha_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                      const shared_ptr<abstract_bipartite_edge> &e,
                                                                      uint32_t alpha)
    {
        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_alpha = get_maximal_b_alpha(B, right_index_map, l, alpha);
        auto tau_alpha = min(left_index_map->at(l)->get_j(alpha),right_index_map->at(r)->get_maximal_j(alpha));
        if(tau_alpha == 0)
        {
            return;
        }

        auto Cl = make_shared<unordered_set<uint32_t>>();
        auto Cr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto sup_l = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto sup_r = make_shared<unordered_map<uint32_t,uint32_t>>();

        if(left_index_map->at(l)->get_j(alpha) == tau_alpha)
        {
            Sl->insert(l);
        }

        if(right_index_map->at(r)->get_maximal_j(alpha) == tau_alpha)
        {
            Sr->insert(r);
        }

        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                if(!B->get_left_vertex(l1)){
                    continue;
                }

                sup_l->insert({l1, 0});
                Tl->insert(l1);

                for(const auto&[r1,e1]:*B->get_left_vertex(l1)->get_edge_map())
                {
                    auto beta_max = right_index_map->at(r1)->get_maximal_j(alpha);
                    if(beta_max >= tau_alpha && ! Cr->count(r1))
                    {
                        ++sup_l->at(l1);
                    }
                }
                if(sup_l->at(l1) < alpha){
                    add_left_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, l1, alpha, tau_alpha);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                sup_r->insert({r1, 0});
                Tr->insert(r1);

                auto r1_vertex = B->get_right_vertex(r1);
                if(!r1_vertex){
                    continue;
                }
                for(const auto& [l1,e1]:*r1_vertex->get_edge_map()){
                    auto beta_max = left_index_map->at(l1)->get_j(alpha);
                    if(beta_max >= tau_alpha && ! Cl->count(l1)){
                        ++sup_r->at(r1);
                    }
                }
                if(sup_r->at(r1) < tau_alpha){
                    add_right_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, r1, alpha, tau_alpha);
                }
            }
        }

        for(const auto&l1:*Cl) {
            left_mutex_map->at(l1)->lock();
            new_left_index_map->at(l1)->remove(alpha, tau_alpha - 1);
            left_mutex_map->at(l1)->unlock();
        }
        left_mutex_map->at(l)->lock();
        new_left_index_map->at(l)->remove(alpha, b_alpha);
        left_mutex_map->at(l)->unlock();

        for(const auto& r1:*Cr){
            uint32_t r1_degree = B->get_right_vertex(r1) ? B->get_right_vertex(r1)->get_degree():0;
            if(r1 == r){
                ++r1_degree;
            }
            right_mutex_map->at(r1)->lock();
            for(uint32_t j = tau_alpha; j<=r1_degree;++j){
                new_right_index_map->at(r1)->remove(j, alpha - 1);
            }
            right_mutex_map->at(r1)->unlock();
        }
        //B->remove_edge(e);
    }

    void edge_bipartite_core_maintenance::update_beta_tau_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                    const shared_ptr<abstract_bipartite_edge>& e,
                                                                    uint32_t beta)
    {
        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_beta = get_maximal_b_beta(B, left_index_map, r, beta);
        auto tau_beta = min(left_index_map->at(l)->get_maximal_i(beta), right_index_map->at(r)->get_i(beta));

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto Cl = make_shared<unordered_set<uint32_t>>();
        auto Cr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto sup_l = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto sup_r = make_shared<unordered_map<uint32_t,uint32_t>>();

        if(right_index_map->at(r)->get_i(beta) == tau_beta)
        {
            Sr->insert(r);
        }

        if(left_index_map->at(l)->get_maximal_i(beta) == tau_beta)
        {
            Sl->insert(l);
        }

        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                sup_l->insert({l1, 0});
                Tl->insert(l1);

                auto l1_vertex = B->get_left_vertex(l1);
                if(!l1_vertex){
                    continue;
                }
                for(const auto&[r1,e1]:*l1_vertex->get_edge_map())
                {
                    auto alpha_max = right_index_map->at(r1)->get_i(beta);
                    if(alpha_max >= tau_beta && ! Cr->count(r1))
                    {
                        ++sup_l->at(l1);
                    }
                }
                if(sup_l->at(l1) < tau_beta)
                {
                    add_left_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, l1, tau_beta, beta);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                sup_r->insert({r1, 0});
                Tr->insert(r1);

                auto r1_vertex = B->get_right_vertex(r1);
                if(!r1_vertex){
                    continue;
                }
                for(const auto& [l1,e1]:*r1_vertex->get_edge_map()){
                    auto alpha_max = left_index_map->at(l1)->get_maximal_i(beta);
                    if(alpha_max >= tau_beta && ! Cl->count(l1)){
                        ++sup_r->at(r1);
                    }
                }
                if(sup_r->at(r1) < beta){
                    add_right_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, r1, tau_beta, beta);
                }
            }
        }

        for(const auto&r1:*Cr)
        {
            if(!new_right_index_map->count(r1)){
                new_right_index_map->insert({r1,make_shared<bipartite_core_right_store_index>()});
            }
            new_right_index_map->at(r1)->insert(beta, tau_beta - 1);
        }

        //B->insert_edge_collection(e);
        for(const auto&l1:*Cl){
            if(!new_left_index_map->count(l1)) {
                new_left_index_map->insert({l1,make_shared<bipartite_core_left_store_index>()});
            }
            uint32_t  l1_degree =  B->get_left_vertex(l1) ? B->get_left_vertex(l1)->get_degree():0;
            if(l1 == l)
            {
                ++l1_degree;
            }
            for(uint32_t i = tau_beta; i <= l1_degree;++i){
                new_left_index_map->at(l1)->remove(i, beta - 1);
            }
        }

        //B->remove_edge(e);

        if(!new_right_index_map->count(r)){
            new_right_index_map->insert({r,make_shared<bipartite_core_right_store_index>()});
        }
        new_right_index_map->at(r)->remove(beta, b_beta);
    }

    void edge_bipartite_core_maintenance::update_beta_tau_beta_core(const shared_ptr<abstract_bipartite_graph> &B,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &left_mutex_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>> &right_mutex_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &right_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_left_store_index>>> &new_left_index_map,
                                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<bipartite_core_right_store_index>>> &new_right_index_map,
                                                                    const shared_ptr<abstract_bipartite_edge> &e,
                                                                    uint32_t beta)
    {
        auto l = e->get_left_vertex_id();
        auto r = e->get_right_vertex_id();
        auto b_beta = get_maximal_b_beta(B, left_index_map, r, beta);
        auto tau_beta = min(left_index_map->at(l)->get_maximal_i(beta), right_index_map->at(r)->get_i(beta));

        auto Tl = make_shared<unordered_set<uint32_t>>();
        auto Tr = make_shared<unordered_set<uint32_t>>();

        auto Cl = make_shared<unordered_set<uint32_t>>();
        auto Cr = make_shared<unordered_set<uint32_t>>();

        auto Sl = make_shared<unordered_set<uint32_t>>();
        auto Sr = make_shared<unordered_set<uint32_t>>();

        auto sup_l = make_shared<unordered_map<uint32_t,uint32_t>>();
        auto sup_r = make_shared<unordered_map<uint32_t,uint32_t>>();

        if(right_index_map->at(r)->get_i(beta) == tau_beta)
        {
            Sr->insert(r);
        }

        if(left_index_map->at(l)->get_maximal_i(beta) == tau_beta)
        {
            Sl->insert(l);
        }

        while (!Sl->empty()||!Sr->empty())
        {
            while(!Sl->empty()){
                auto l1 = *Sl->begin();
                Sl->erase(l1);

                sup_l->insert({l1, 0});
                Tl->insert(l1);

                auto l1_vertex = B->get_left_vertex(l1);
                if(!l1_vertex){
                    continue;
                }
                for(const auto&[r1,e1]:*l1_vertex->get_edge_map())
                {
                    auto alpha_max = right_index_map->at(r1)->get_i(beta);
                    if(alpha_max >= tau_beta && ! Cr->count(r1))
                    {
                        ++sup_l->at(l1);
                    }
                }
                if(sup_l->at(l1) < tau_beta)
                {
                    add_left_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, l1, tau_beta, beta);
                }
            }
            while(!Sr->empty())
            {
                auto r1 = *Sr->begin();
                Sr->erase(r1);

                sup_r->insert({r1, 0});
                Tr->insert(r1);

                auto r1_vertex = B->get_right_vertex(r1);
                if(!r1_vertex){
                    continue;
                }
                for(const auto& [l1,e1]:*r1_vertex->get_edge_map()){
                    auto alpha_max = left_index_map->at(l1)->get_maximal_i(beta);
                    if(alpha_max >= tau_beta && ! Cl->count(l1)){
                        ++sup_r->at(r1);
                    }
                }
                if(sup_r->at(r1) < beta){
                    add_right_candidates(B, left_index_map, right_index_map, Tl, Tr, Cl, Cr, Sl, Sr, sup_l, sup_r, r1, tau_beta, beta);
                }
            }
        }

        for(const auto&r1:*Cr)
        {
            right_mutex_map->at(r1)->lock();
            new_right_index_map->at(r1)->insert(beta, tau_beta - 1);
            right_mutex_map->at(r1)->unlock();
        }
        right_mutex_map->at(r)->lock();
        new_right_index_map->at(r)->remove(beta, b_beta);
        right_mutex_map->at(r)->unlock();


        for(const auto&l1:*Cl){
            uint32_t  l1_degree = B->get_left_vertex(l1)->get_degree();
            if(l1 == l)
            {
                ++l1_degree;
            }
            left_mutex_map->at(l1)->lock();
            for(uint32_t i = tau_beta; i <= l1_degree;++i){
                new_left_index_map->at(l1)->remove(i, beta - 1);
            }
            left_mutex_map->at(l1)->unlock();
        }
    }
}