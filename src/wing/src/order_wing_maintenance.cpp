
#include "wing/order_wing_maintenance.h"

namespace scnu{

    void order_wing_maintenance::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &C_k,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& s,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &evicted_edge_set,
                                                          uint32_t k){
        for(const auto& e1:*C_k)
        {
            s->insert({e1,0});
            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            for(const auto&[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()){
                if(r2 == r1 || !(C_k->count(e2) || edge_wing_map->at(e2) >= k)){
                    continue;
                }

                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || !(C_k->count(e3) || edge_wing_map->at(e3)>= k)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || !(C_k->count(e4) || edge_wing_map->at(e4) >= k)) {
                        continue;
                    }

                    ++s->at(e1);
                }
            }

            if(s->at(e1) <= k){
                evicted_edge_set->insert(e1);
            }
        }
    }


    void order_wing_maintenance::init(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &ext){

        for(const auto&[e, wing_number]:*edge_wing_map)
        {
            ext->insert({e,0});
        }
    }

    void order_wing_maintenance::insert(const shared_ptr<abstract_bipartite_graph>& G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_wing_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& ts,
                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<double,shared_ptr<abstract_bipartite_edge>>>>>& wing_order_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& rem,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& ext) {
        for(const auto &e:*edge_set)
        {
            G->insert_edge(e);
            rem->insert({e,0});
            ext->insert({e,0});
            ts->insert({e,0});
            edge_wing_map->insert({e,0});
        }

        /**
         * @brief store edges need tp update their ts
         */
        auto N_k = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(edge_set);

        auto evicted_edge_set =  make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for(uint32_t  k = 0;!N_k->empty(); ++k)
        {
            if(!wing_order_map->count(k)){
                wing_order_map->insert({k,make_shared<extend_list<double,shared_ptr<abstract_bipartite_edge>>>()});
            }

            auto current_order_list = wing_order_map->at(k);

            auto C_k = container_copy::to_unordered_set<shared_ptr<abstract_bipartite_edge>>(N_k);

            auto s = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>();
            auto e_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();
            update_order_list(G, evicted_edge_set, C_k, rem, s, edge_wing_map, current_order_list, e_vector, k);

            auto B = make_shared<map<double,shared_ptr<extend_node<double,shared_ptr<abstract_bipartite_edge>>>>>();
            update_edge_support(G,C_k,edge_wing_map,wing_order_map,ext,B,k);

            auto p = current_order_list->get_head();

            for(auto iter = e_vector->rbegin(); iter!=e_vector->rend(); ++iter){
                auto e_node = make_shared<extend_node<double,shared_ptr<abstract_bipartite_edge>>>(0, *iter);
                current_order_list->left_insert(e_node);
            }

            while(p){
                auto p_next = p->get_next();
                auto e_current = p->get_value();

                B->erase(p->get_key());

                if(ext->at(e_current) == 0){
                    p_next = B->empty() ? shared_ptr<extend_node<double,shared_ptr<abstract_bipartite_edge>>>(): B->begin()->second;
                }else if(rem->at(e_current)+ext->at(e_current) > k){

                    C_k->insert(e_current);

                    s->insert({e_current,ext->at(e_current)+rem->at(e_current)});

                    ext->at(e_current) = 0;

                    update_edge_support(G,C_k,edge_wing_map,wing_order_map,ext,B,e_current,k);


                    current_order_list->remove(e_current);
                }else
                {
                    rem->at(e_current) = ext->at(e_current)+rem->at(e_current);
                    ext->at(e_current) = 0;

                    auto l1 = e_current->get_left_vertex_id();
                    auto r1 = e_current->get_right_vertex_id();

                    for(const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()){
                        if(r2 == r1 || !(C_k->count(e2) || edge_wing_map->at(e2) > k || current_order_list->count_value(e2) && test_order(edge_wing_map, wing_order_map, e_current,e2))) {
                            continue;
                        }
                        for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                            if (l2 == l1 ||!(C_k->count(e3) || edge_wing_map->at(e3) > k || current_order_list->count_value(e3) && test_order(edge_wing_map, wing_order_map, e_current,e3))) {
                                continue;
                            }
                            auto e4 = G->get_edge(l2, r2);
                            if (!e4 || !(C_k->count(e4) || edge_wing_map->at(e4) > k || current_order_list->count_value(e4) && test_order(edge_wing_map, wing_order_map, e_current,e4))) {
                                continue;
                            }

                            if(!C_k->count(e2) && !C_k->count(e3) && !C_k->count(e4)){
                                continue;
                            }

                            if (C_k->count(e2)) {
                                --s->at(e2);
                                if (s->at(e2) <= k) {
                                    evicted_edge_set->insert(e2);
                                }
                            }else if(ext->at(e2) > 0){
                                --ext->at(e2);
                                if (ext->at(e2) == 0) {
                                    auto e2_node = current_order_list->find(e2);
                                    B->erase(e2_node->get_key());
                                }
                            }


                            if (C_k->count(e3)) {
                                --s->at(e3);
                                if (s->at(e3) <= k) {
                                    evicted_edge_set->insert(e3);
                                }
                            }else if(ext->at(e3) > 0){
                                --ext->at(e3);
                                if (ext->at(e3) == 0) {
                                    auto e3_node = current_order_list->find(e3);
                                    B->erase(e3_node->get_key());
                                }
                            }

                            if (C_k->count(e4)) {
                                --s->at(e4);
                                if (s->at(e4) <= k) {
                                    evicted_edge_set->insert(e4);
                                }
                            }else if(ext->at(e4) > 0){
                                --ext->at(e4);
                                if (ext->at(e4) == 0) {
                                    auto e4_node = current_order_list->find(e4);
                                    B->erase(e4_node->get_key());
                                }
                            }
                        }
                    }
                    remove_unsatisfied_edges(G,evicted_edge_set,C_k,edge_wing_map,wing_order_map,s,rem,ext,B, p_next,e_current,k);
                }
                p = p_next;
            }

            if(!current_order_list->empty()){
                current_order_list->reset_order();
            }

            auto ts_update = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            for(const auto&e:*N_k){
                if(!C_k->count(e)){
                    edge_wing_map->at(e) = k;
                    ts_update->insert(e);
                    ts->at(e) = 0;
                }
            }
            increase_ts(G, ts_update, edge_wing_map, C_k, ts, k);
            N_k = C_k;
        }
    }

    uint32_t order_wing_maintenance::k_new_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                       const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                       const shared_ptr<abstract_bipartite_edge> &e)
    {
        const auto&e1 = e;
        auto result_map = make_shared<map<uint32_t, uint32_t>>();

        auto l1 = e1->get_left_vertex_id();
        auto r1 = e1->get_right_vertex_id();

        for(const auto&[l2,e2]:*G->get_right_vertex(r1)->get_edge_map()){
            if(l2 == l1){
                continue;
            }

            for(const auto&[r2,e3]:*G->get_left_vertex(l1)->get_edge_map()){
                if(r2 == r1){
                    continue;
                }
                auto e4 = G->get_edge(l2,r2);
                if(!e4){
                    continue;
                }

                auto min_k = min({edge_wing_map->at(e2),edge_wing_map->at(e3),edge_wing_map->at(e4)});
                if(!result_map->count(min_k)){
                    result_map->insert({min_k, 0});
                }
                ++result_map->at(min_k);
            }
        }

        uint32_t  k = 0;
        if(!result_map->empty()){
            uint32_t sum = 0;
            k = result_map->rbegin()->first;
            while(k >= 0)
            {
                if(result_map->count(k)){
                    sum += result_map->at(k);
                }
                if(sum >= k){
                    break;
                }
                --k;
            }
        }
        return k;
    }

    void order_wing_maintenance::remove(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &ts,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem) {

        auto Q = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        update_rem_and_ts(G,edge_set,edge_wing_map,wing_order_map,rem,ts,Q);

        for(const auto&e1:*edge_set){
            G->remove_edge(e1);
            auto k = edge_wing_map->at(e1);
            wing_order_map->at(k)->remove(e1);
            edge_wing_map->erase(e1);
        }

        while (!Q->empty()) {
            auto e1 = *Q->begin();
            Q->erase(e1);

            auto k = edge_wing_map->at(e1);
            auto k_new = k_new_computation(G, edge_wing_map, e1);

            rem->at(e1) = 0;

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || edge_wing_map->at(e2) <= k_new) {
                    continue;
                }
                for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || edge_wing_map->at(e3) <= k_new) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || edge_wing_map->at(e4) <= k_new) {
                        continue;
                    }

                    ++rem->at(e1);

                    auto min_e = e1;

                    if (test_order(edge_wing_map, wing_order_map, e2, min_e)) {
                        min_e = e2;
                    }

                    if (test_order(edge_wing_map, wing_order_map, e3, min_e)) {
                        min_e = e3;
                    }

                    if (test_order(edge_wing_map, wing_order_map, e4, min_e)) {
                        min_e = e4;
                    }

                    if (min_e == e2) {
                        --rem->at(e2);
                    }

                    if (min_e == e3) {
                        --rem->at(e3);
                    }

                    if (min_e == e4) {
                        --rem->at(e4);
                    }

                    auto min_k = min(edge_wing_map->at(e2), edge_wing_map->at(e3));
                    min_k = min(min_k, edge_wing_map->at(e4));
                    if (min_k <= edge_wing_map->at(e1)) {
                        if (edge_wing_map->at(e2) == min_k) {
                            --ts->at(e2);
                            if (ts->at(e2) < edge_wing_map->at(e2) && !Q->count(e2)) {
                                Q->insert(e2);
                            }
                        }

                        if (edge_wing_map->at(e3) == min_k) {
                            --ts->at(e3);
                            if (ts->at(e3) < edge_wing_map->at(e3) && !Q->count(e3)) {
                                Q->insert(e3);
                            }
                        }

                        if (edge_wing_map->at(e4) == min_k) {
                            --ts->at(e4);
                            if (ts->at(e4) < edge_wing_map->at(e4) && !Q->count(e4)) {
                                Q->insert(e4);
                            }
                        }
                    }
                }
            }

            wing_order_map->at(k)->remove(e1);
            if (!wing_order_map->count(k_new)) {
                wing_order_map->insert(
                        {k_new, make_shared<extend_list<double, shared_ptr<abstract_bipartite_edge>>>()});
            }
            wing_order_map->at(k_new)->push_back(e1);
            edge_wing_map->at(e1) = k_new;


            /**
             * @brief recompute ts
             */
            ts->at(e1) = 0;
            for (const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || edge_wing_map->at(e2) < k_new) {
                    continue;
                }
                for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || edge_wing_map->at(e3) < k_new) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || edge_wing_map->at(e4) < k_new) {
                        continue;
                    }
                    ++ts->at(e1);
                }
            }
        }
    }


    void order_wing_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_bipartite_graph> &G,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& C_k,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &s,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &ext,
                                                          const shared_ptr<map<double,shared_ptr<extend_node<double,shared_ptr<abstract_bipartite_edge>>>>>& B,
                                                          const shared_ptr<extend_node<double, shared_ptr<abstract_bipartite_edge>>>& e_pivot_node,
                                                          const shared_ptr<abstract_bipartite_edge>& e_current,
                                                          uint32_t k) {
        auto current_order_list = wing_order_map->at(k);
        auto e_vector = make_shared<vector<shared_ptr<abstract_bipartite_edge>>>();

        while(!current_edge_set->empty())
        {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            for(const auto &e1:*current_edge_set){
                C_k->erase(e1);
                s->erase(e1);

                rem->at(e1) = 0;

                e_vector->push_back(e1);

                auto l1 = e1->get_left_vertex_id();
                auto r1 = e1->get_right_vertex_id();
                for(const auto&[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                    if (r2 == r1 || !(C_k->count(e2) || edge_wing_map->at(e2) > k || current_order_list->count_value(e2) && test_order(edge_wing_map, wing_order_map,  e_current ,e2))) {
                        continue;
                    }
                    for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                        if (l2 == l1 || !(C_k->count(e3) || edge_wing_map->at(e3) > k || current_order_list->count_value(e3) && test_order(edge_wing_map, wing_order_map,  e_current ,e3))) {
                            continue;
                        }

                        auto e4 = G->get_edge(l2, r2);
                        if (!e4 || !(C_k->count(e4) || edge_wing_map->at(e4) > k || current_order_list->count_value(e4) && test_order(edge_wing_map, wing_order_map,  e_current ,e4))) {
                            continue;
                        }

                        ++rem->at(e1);

                        if (C_k->count(e2)) {
                            --s->at(e2);
                            if (s->at(e2) <= k  && !current_edge_set->count(e2)) {
                                next_edge_set->insert(e2);
                            }
                        } else if (ext->at(e2) > 0) {
                            --ext->at(e2);
                            if (ext->at(e2) == 0) {
                                auto e2_node = current_order_list->find(e2);
                                B->erase(e2_node->get_key());
                            }
                        }

                        if (C_k->count(e3)) {
                            --s->at(e3);
                            if (s->at(e3) <= k  && !current_edge_set->count(e3)) {
                                next_edge_set->insert(e3);
                            }
                        } else if (ext->at(e3) > 0) {
                            --ext->at(e3);
                            if (ext->at(e3) == 0) {
                                auto e3_node = current_order_list->find(e3);
                                B->erase(e3_node->get_key());
                            }
                        }

                        if (C_k->count(e4)) {
                            --s->at(e4);
                            if (s->at(e4) <= k  && !current_edge_set->count(e4)) {
                                next_edge_set->insert(e4);
                            }
                        } else if (ext->at(e4) > 0) {
                            --ext->at(e4);
                            if (ext->at(e4) == 0) {
                                auto e4_node = current_order_list->find(e4);
                                B->erase(e4_node->get_key());
                            }
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*e_vector){
            auto e_node = make_shared<extend_node<double,shared_ptr<abstract_bipartite_edge>>>(0, e);
            if(!e_pivot_node){
                current_order_list->right_insert(e_node);
            }else{
                current_order_list->insert_before(e_node, e_pivot_node);
            }
        }
    }

    bool order_wing_maintenance::test_order(
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
            const shared_ptr<abstract_bipartite_edge> &e1,
            const shared_ptr<abstract_bipartite_edge> &e2) {
        auto k1 = edge_wing_map->at(e1);
        auto k2 = edge_wing_map->at(e2);
        if(k1<k2)
        {
            return  true;
        }
        else if(k1>k2){
            return false;
        }else
        {
            auto order_list = wing_order_map->at(k1);
            return order_list->find_key(e1).value() < order_list->find_key(e2).value();
        }
    }



    /**
     * @details update the support affected by e
     * @param G
     * @param C_k
     * @param edge_wing_map
     * @param wing_order_map
     * @param ext
     * @param B
     * @param e
     * @param k
     */
    void order_wing_maintenance::update_edge_support(const shared_ptr<abstract_bipartite_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& C_k,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& ext,
                                                     const shared_ptr<map<double,shared_ptr<extend_node<double,shared_ptr<abstract_bipartite_edge>>>>>& B,
                                                     const shared_ptr<abstract_bipartite_edge> &e,
                                                     uint32_t k) {
        auto current_order_list = wing_order_map->at(k);

        const auto &e1 = e;
        auto l1 = e1->get_left_vertex_id();
        auto r1 = e1->get_right_vertex_id();

        for(const auto &[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()) {
            if (r2 == r1 || C_k->count(e2) || test_order(edge_wing_map, wing_order_map, e2, e1)) {
                continue;
            }
            for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                if (l2 == l1 || C_k->count(e3) || test_order(edge_wing_map, wing_order_map, e3, e1)) {
                    continue;
                }

                auto e4 = G->get_edge(l2, r2);
                if (!e4 || C_k->count(e4) || test_order(edge_wing_map, wing_order_map, e4, e1)) {
                    continue;
                }


                if (current_order_list->count_value(e2)) {
                    ++ext->at(e2);
                    auto e2_node = current_order_list->find(e2);
                    B->insert({e2_node->get_key(), e2_node});
                }

                if (current_order_list->count_value(e3)) {
                    ++ext->at(e3);
                    auto e3_node = current_order_list->find(e3);
                    B->insert({e3_node->get_key(), e3_node});
                }


                if (current_order_list->count_value(e4)) {
                    ++ext->at(e4);
                    auto e4_node = current_order_list->find(e4);
                    B->insert({e4_node->get_key(), e4_node});
                }
            }
        }
    }

    /**
     * @details update the support affected by edges in C_k
     * @param G
     * @param C_k
     * @param edge_wing_map
     * @param wing_order_map
     * @param ext
     * @param B
     * @param e
     * @param k
     */
    void order_wing_maintenance::update_edge_support(const shared_ptr<abstract_bipartite_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &C_k,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& ext,
                                                     const shared_ptr<map<double,shared_ptr<extend_node<double,shared_ptr<abstract_bipartite_edge>>>>>& B,
                                                     uint32_t k) {
        auto current_order_list = wing_order_map->at(k);

        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        for(const auto&e1:*C_k){
            visited_edge_set->insert(e1);

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || visited_edge_set->count(e2) || !(C_k->count(e2) || edge_wing_map->at(e2) >= k)) {
                    continue;
                }

                for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || visited_edge_set->count(e3) || !(C_k->count(e3) || edge_wing_map->at(e3) >= k)) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || visited_edge_set->count(e4) || !(C_k->count(e4) || edge_wing_map->at(e4) >= k)) {
                        continue;
                    }

                    if (current_order_list->count_value(e2)) {
                        ++ext->at(e2);
                        auto e2_node = current_order_list->find(e2);
                        B->insert({e2_node->get_key(), e2_node});
                    }


                    if (current_order_list->count_value(e3)) {
                        ++ext->at(e3);
                        auto e3_node = current_order_list->find(e3);
                        B->insert({e3_node->get_key(), e3_node});
                    }


                    if (current_order_list->count_value(e4)) {
                        ++ext->at(e4);
                        auto e4_node = current_order_list->find(e4);
                        B->insert({e4_node->get_key(), e4_node});
                    }
                }
            }
        }
    }

    void order_wing_maintenance::decrease_rem(const shared_ptr<abstract_bipartite_graph> &G,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                              const shared_ptr<abstract_bipartite_edge> &e) {


        const auto &e1 = e;

        rem->at(e1) = 0;
        auto l1 = e1->get_left_vertex_id();
        auto r1 = e1->get_right_vertex_id();

        if(!G->get_left_vertex(l1) || !G->get_right_vertex(r1)){
            return;
        }

        for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
            if (r2 == r1) {
                continue;
            }
            for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                if (l2 == l1) {
                    continue;
                }
                auto e4 = G->get_edge(l2, r2);
                if (!e4) {
                    continue;
                }

                auto min_e = e1;

                if(test_order(edge_wing_map, wing_order_map, e2, min_e))
                {
                    min_e = e2;
                }

                if(test_order(edge_wing_map, wing_order_map, e3, min_e)){
                    min_e = e3;
                }

                if(test_order(edge_wing_map, wing_order_map, e4, min_e)){
                    min_e = e4;
                }

                if (e2 == min_e) {
                    --rem->at(e2);
                }

                if (e3 == min_e) {
                    --rem->at(e3);
                }

                if (e4 == min_e) {
                    --rem->at(e4);
                }
            }
        }
    }

    void order_wing_maintenance::increase_rem(const shared_ptr<abstract_bipartite_graph> &G,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                              const shared_ptr<abstract_bipartite_edge> &e) {


        const auto &e1 = e;

        rem->at(e1) = 0;
        auto l1 = e1->get_left_vertex_id();
        auto r1 = e1->get_right_vertex_id();

        if(!G->get_left_vertex(l1) || !G->get_right_vertex(r1)){
            return;
        }

        for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
            if (r2 == r1) {
                continue;
            }
            for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                if (l2 == l1) {
                    continue;
                }
                auto e4 = G->get_edge(l2, r2);
                if (!e4) {
                    continue;
                }

                auto min_e = e1;
                if(test_order(edge_wing_map, wing_order_map, e2, min_e)){
                    min_e = e2;
                }

                if(test_order(edge_wing_map, wing_order_map, e3, min_e)){
                    min_e = e3;
                }

                if(test_order(edge_wing_map, wing_order_map, e4, min_e)){
                    min_e = e4;
                }

                if(e1 == min_e){
                    ++rem->at(e1);
                }

                if (e2 == min_e) {
                    ++rem->at(e2);
                }

                if (e3 == min_e) {
                    ++rem->at(e3);
                }

                if (e4 == min_e) {
                    ++rem->at(e4);
                }
            }
        }
    }

    void order_wing_maintenance::update_rem_and_ts(const shared_ptr<abstract_bipartite_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> & edge_set,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>>>> &wing_order_map,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &ts,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &Q) {
        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        for(const auto &e1:*edge_set){
            visited_edge_set->insert(e1);

            ts->erase(e1);
            rem->erase(e1);

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || visited_edge_set->count(e2)) {
                    continue;
                }
                for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || visited_edge_set->count(e3)) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || visited_edge_set->count(e4)) {
                        continue;
                    }


                    if (edge_wing_map->at(e2) > edge_wing_map->at(e1) &&
                        edge_wing_map->at(e3) > edge_wing_map->at(e1) &&
                        edge_wing_map->at(e4) > edge_wing_map->at(e1)) {
                        continue;
                    }

                    auto min_k = min(edge_wing_map->at(e2), edge_wing_map->at(e3));
                    min_k = min(min_k, edge_wing_map->at(e4));

                    if (edge_wing_map->at(e2) == min_k) {
                        --ts->at(e2);
                        if (!edge_set->count(e2) && ts->at(e2) < edge_wing_map->at(e2)) {
                            Q->insert(e2);
                        }
                    }

                    if (edge_wing_map->at(e3) == min_k) {
                        --ts->at(e3);
                        if (!edge_set->count(e3) && ts->at(e3) < edge_wing_map->at(e3)) {
                            Q->insert(e3);
                        }
                    }

                    if (edge_wing_map->at(e4) == min_k) {
                        --ts->at(e4);
                        if (!edge_set->count(e4) && ts->at(e4) < edge_wing_map->at(e4)) {
                            Q->insert(e4);
                        }
                    }
                    auto min_e = e1;

                    if (test_order(edge_wing_map, wing_order_map, e2, min_e)) {
                        min_e = e2;
                    }

                    if (test_order(edge_wing_map, wing_order_map, e3, min_e)) {
                        min_e = e3;
                    }

                    if (test_order(edge_wing_map, wing_order_map, e4, min_e)) {
                        min_e = e4;
                    }

                    if (min_e == e2) {
                        --rem->at(e2);
                    }

                    if (min_e == e3) {
                        --rem->at(e3);
                    }

                    if (min_e == e4) {
                        --rem->at(e4);
                    }
                }
            }
        }
    }



    void order_wing_maintenance::update_order_list(const shared_ptr<abstract_bipartite_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &current_edge_set,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &C_k,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &rem,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &s,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                   const shared_ptr<extend_list<double, shared_ptr<abstract_bipartite_edge>>> &current_order_list,
                                                   const shared_ptr<vector<shared_ptr<abstract_bipartite_edge>>> &e_vector,
                                                   uint32_t k){
        for(const auto& e1:*C_k)
        {
            s->insert({e1,0});
            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            for(const auto&[r2,e2]:*G->get_left_vertex(l1)->get_edge_map()){
                if(r2 == r1 || !C_k->count(e2) && edge_wing_map->at(e2) < k){
                    continue;
                }

                for (const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || !C_k->count(e3) && edge_wing_map->at(e3) < k) {
                        continue;
                    }

                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || !C_k->count(e4) && edge_wing_map->at(e4) < k) {
                        continue;
                    }

                    ++s->at(e1);
                }
            }

            if(s->at(e1) <= k){
                current_edge_set->insert(e1);
            }
        }

        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

            for(const auto &e1:*current_edge_set){
                C_k->erase(e1);
                rem->at(e1) = 0;

                e_vector->push_back(e1);


                auto l1 = e1->get_left_vertex_id();
                auto r1 = e1->get_right_vertex_id();

                for (const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                    if (r2 == r1 || !(C_k->count(e2) || edge_wing_map->at(e2) >= k)) {
                        continue;
                    }

                    for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                        if (l2 == l1 || !(C_k->count(e3) || edge_wing_map->at(e3) >= k)) {
                            continue;
                        }

                        auto e4 = G->get_edge(l2, r2);
                        if (!e4 || !(C_k->count(e4) || edge_wing_map->at(e4) >= k)) {
                            continue;
                        }

                        ++rem->at(e1);

                        if (C_k->count(e2)) {
                            --s->at(e2);
                            if (s->at(e2) <= k && !current_edge_set->count(e2)) {
                                next_edge_set->insert(e2);
                            }
                        }

                        if (C_k->count(e3)) {
                            --s->at(e3);
                            if (s->at(e3) <= k && !current_edge_set->count(e3)) {
                                next_edge_set->insert(e3);
                            }
                        }

                        if (C_k->count(e4)) {
                            --s->at(e4);
                            if (s->at(e4) <= k && !current_edge_set->count(e4)) {
                                next_edge_set->insert(e4);
                            }
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);


        }


    }

    void order_wing_maintenance::increase_ts(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &C_k,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &ts,
                                             uint32_t k) {
        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
        /*while(!edge_set->empty()){
            auto e1 = *edge_set->begin();
            edge_set->erase(e1);

            if(visited_edge_set->count(e1)){
                continue;
            }
            visited_edge_set->insert(e1);

            ts->at(e1) = 0;

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            for(const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1) {
                    continue;
                }
                for (const auto&[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4) {
                        continue;
                    }

                    if(edge_wing_map->at(e2) >= edge_wing_map->at(e1)
                       && edge_wing_map->at(e3) >= edge_wing_map->at(e1)
                       && edge_wing_map->at(e4) >= edge_wing_map->at(e1)){
                        ++ts->at(e1);
                    }

                    if(edge_wing_map->at(e2) <=  edge_wing_map->at(e1)){
                        edge_set->insert(e2);
                    }

                    if(edge_wing_map->at(e3) <=  edge_wing_map->at(e1)){
                        edge_set->insert(e3);
                    }

                    if(edge_wing_map->at(e4) <=  edge_wing_map->at(e1)){
                        edge_set->insert(e4);
                    }
                }
            }
        }*/

        for (const auto &e1: *edge_set) {
            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();

            visited_edge_set->insert(e1);



            for(const auto&[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()) {
                if (r2 == r1 || visited_edge_set->count(e2) || !(C_k->count(e2) || edge_wing_map->at(e2) >= k)) {
                    continue;
                }
                for (const auto &[l2, e3]: *G->get_right_vertex(r1)->get_edge_map()) {
                    if (l2 == l1 || visited_edge_set->count(e3) || !(C_k->count(e3) || edge_wing_map->at(e3) >= k)) {
                        continue;
                    }
                    auto e4 = G->get_edge(l2, r2);
                    if (!e4 || visited_edge_set->count(e4) || !(C_k->count(e4) || edge_wing_map->at(e4) >= k)) {
                        continue;
                    }

                    ++ts->at(e1);

                    if (edge_wing_map->at(e2) == k) {
                        ++ts->at(e2);
                    }

                    if (edge_wing_map->at(e3) == k) {
                        ++ts->at(e3);
                    }

                    if (edge_wing_map->at(e4) == k) {
                        ++ts->at(e4);
                    }
                }
            }
        }
    }
}