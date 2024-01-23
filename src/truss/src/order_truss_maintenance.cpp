
#include "truss/order_truss_maintenance.h"

namespace scnu{

    void order_truss_maintenance::edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& s,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                                          uint32_t k){

        auto affected_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto& e1:*C_k)
        {
            s->insert({e1,0});
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                swap(u,v);
            }

            for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v || !(C_k->count(e2) || edge_truss_map->at(e2) >= k)){
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || !(C_k->count(e3) || edge_truss_map->at(e3) >= k)) {
                    continue;
                }

                ++s->at(e1);

                if(edge_truss_map->at(e2) == k){
                    affected_set->insert(e2);
                }

                if(edge_truss_map->at(e3) == k){
                    affected_set->insert(e3);
                }
            }

            if (s->at(e1) <= k - 2) {
                evicted_edge_set->insert(e1);
            }
        }

        for(const auto &e1:*affected_set){
            s->insert({e1,0});
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                swap(u,v);
            }

            for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v || !(C_k->count(e2) || edge_truss_map->at(e2) >= k)){
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || !(C_k->count(e3) || edge_truss_map->at(e3) >= k)) {
                    continue;
                }

                ++s->at(e1);
            }

            if (s->at(e1) <= k - 2) {
                evicted_edge_set->insert(e1);
            }
        }
    }


    uint32_t order_truss_maintenance::compute_edge_support(const shared_ptr<abstract_graph> &G,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                           const shared_ptr<abstract_edge> &e1,
                                                           const shared_ptr<abstract_edge> &e_pivot,
                                                           uint32_t k){
        uint32_t support = 0;

        auto order_list = truss_order_map->at(k);

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
            swap(u, v);
        }

        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
            if (w == v || !(C_k->count(e2) || edge_truss_map->at(e2) > k ||
                            order_list->count_value(e2) && test_order(order_list, e_pivot, e2))) {
                continue;
            }
            auto e3 = G->get_edge(v, w);
            if (!e3 || !(C_k->count(e3) || edge_truss_map->at(e3) > k ||
                         order_list->count_value(e3) && test_order(order_list, e_pivot, e3))) {
                continue;
            }

            ++support;
        }
        return support;
    }


    void
    order_truss_maintenance::init(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                  const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext) {
        for (const auto &[e, truss_number]: *edge_truss_map) {
            ext->insert({e, 0});
        }
    }

    void order_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& edge_truss_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& ts,
                                        const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<int,shared_ptr<abstract_edge>>>>>& truss_order_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& rem,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& ext) {
        for(const auto &e:*edge_set)
        {
            G->insert_edge(e);
            rem->insert({e,0});
            ext->insert({e,0});
            ts->insert({e,0});
            edge_truss_map->insert({e, 2});
        }

        /**
         * @brief store edges need tp update their ts
         */
        auto N_k = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(edge_set);

        auto evicted_edge_set =  make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        for(uint32_t  k = 2;!N_k->empty(); ++k)
        {
            if(!truss_order_map->count(k)){
                truss_order_map->insert({k, make_shared<extend_list<int,shared_ptr<abstract_edge>>>()});
            }

            auto current_order_list = truss_order_map->at(k);

//            for (auto p = current_order_list->get_head(); p; p = p->get_next()) {
//                auto e = p->get_value();
//                printf("(%u, %u) ", e->get_source_vertex_id(), e->get_destination_vertex_id());
//            }
//            printf("\n");

            auto C_k = container_copy::to_unordered_set<shared_ptr<abstract_edge>>(N_k);

            auto s = make_shared<unordered_map<shared_ptr<abstract_edge>,uint32_t>>();
            auto e_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
            update_order_list(G, evicted_edge_set, C_k, rem, s, edge_truss_map, current_order_list,e_vector, k);

            auto B = make_shared<map<int,shared_ptr<extend_node<int,shared_ptr<abstract_edge>>>>>();
            update_edge_support(G, C_k, edge_truss_map, truss_order_map, ext, B, k);

            auto p = current_order_list->get_head();

            for(auto iter = e_vector->rbegin(); iter!=e_vector->rend(); ++iter){
                auto e_node = make_shared<extend_node<int, shared_ptr<abstract_edge>>>(0, *iter);
                current_order_list->left_insert(e_node);
            }

            while(p){
                auto p_next = p->get_next();
                auto e_current = p->get_value();

                B->erase(p->get_key());

                if(ext->at(e_current) == 0){
                    p_next = B->empty() ? shared_ptr<extend_node<int,shared_ptr<abstract_edge>>>(): B->begin()->second;
                }else if(ext->at(e_current) + rem->at(e_current) > k - 2){
                    /**
                     * @brief insert e_current into C_k
                     */
                    C_k->insert(e_current);

                    s->insert({e_current, ext->at(e_current) + rem->at(e_current)});

                    ext->at(e_current) = 0;

                    update_edge_support(G, C_k, edge_truss_map, truss_order_map, ext, B, e_current, k);

                    /**
                     * @brief remove e_current from wing order
                     */
                    current_order_list->remove(e_current);
                }else
                {
                    rem->at(e_current) = ext->at(e_current) + rem->at(e_current);
                    ext->at(e_current) = 0;

                    auto u = e_current->get_source_vertex_id();
                    auto v = e_current->get_destination_vertex_id();

                    if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                        swap(u,v);
                    }

                    for(const auto&[w, e2]:*G->get_vertex(u)->get_edge_map()){
                        if(w == v || !(C_k->count(e2) || edge_truss_map->at(e2) > k ||current_order_list->count_value(e2) &&
                                                                                              test_order(edge_truss_map, truss_order_map, e_current, e2))){
                            continue;
                        }
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || !(C_k->count(e3) || edge_truss_map->at(e3) > k||current_order_list->count_value(e3) &&
                                                                                                      test_order(edge_truss_map, truss_order_map, e_current, e3))) {
                            continue;
                        }

                        if (!C_k->count(e2) && !C_k->count(e3)) {
                            continue;
                        }

                        if (C_k->count(e2)) {
                            --s->at(e2);
                            if (s->at(e2) <= k - 2) {
                                evicted_edge_set->insert(e2);
                            }
                        }else if (ext->at(e2) > 0) {
                            --ext->at(e2);
                            if (ext->at(e2) == 0) {
                                auto e2_node = current_order_list->find(e2);
                                B->erase(e2_node->get_key());
                            }
                        }

                        if (C_k->count(e3)) {
                            --s->at(e3);
                            if (s->at(e3) <= k - 2) {
                                evicted_edge_set->insert(e3);
                            }
                        }else if (ext->at(e3) > 0) {
                            --ext->at(e3);
                            if (ext->at(e3) == 0) {
                                auto e3_node = current_order_list->find(e3);
                                B->erase(e3_node->get_key());
                            }
                        }
                    }

                    remove_unsatisfied_edges(G, evicted_edge_set, C_k, edge_truss_map, truss_order_map, s, rem, ext, B,  p_next, e_current, k);
               }
                p = p_next;
            }

            if(!current_order_list->empty()){
                current_order_list->reset_order();
            }

            auto ts_update = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto&e:*N_k){
                if(!C_k->count(e)){
                    edge_truss_map->at(e) = k;
                    ts_update->insert(e);
                    ts->at(e) = 0;
                }
            }
            increase_ts(G, ts_update, edge_truss_map, C_k, ts, k);
            N_k = C_k;
        }
    }


    uint32_t order_truss_maintenance::k_new_computation(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                        const shared_ptr<abstract_edge> &e) {
        const auto &e1 = e;
        auto result_map = make_shared<map<uint32_t, uint32_t>>();

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
            swap(u,v);
        }

        for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
            if(w == v){
                continue;
            }

            auto e3 = G->get_edge(v, w);
            if (!e3) {
                continue;
            }

            auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
            if(!result_map->count(min_k)){
                result_map->insert({min_k, 0});
            }
            ++result_map->at(min_k);
        }


        uint32_t k = 2;
        if(!result_map->empty()){
            uint32_t sum = 0;
            k = result_map->rbegin()->first;
            while(k >= 2)
            {
                if(result_map->count(k)){
                    sum += result_map->at(k);
                }
                if(sum >= k - 2){
                    break;
                }
                --k;
            }
        }
        return k;
    }

    void order_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ts,
                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem) {

        auto Q = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        update_rem_and_ts(G, edge_set, edge_truss_map, truss_order_map, rem, ts, Q);

        for (const auto &e1: *edge_set) {
            G->remove_edge(e1);
            auto k = edge_truss_map->at(e1);
            truss_order_map->at(k)->remove(e1);
            edge_truss_map->erase(e1);
        }

        auto k_map = make_shared<map<uint32_t, shared_ptr<vector<shared_ptr<scnu::abstract_edge>>>>>();

        while (!Q->empty()) {
            auto e1 = *Q->begin();
            Q->erase(e1);

            auto k = edge_truss_map->at(e1);
            auto k_new = k_new_computation(G, edge_truss_map, e1);

            rem->at(e1) = 0;

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();


            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || edge_truss_map->at(e2) <= k_new) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_truss_map->at(e3) <= k_new) {
                    continue;
                }

                ++rem->at(e1);

                auto min_e = e1;

                if (test_order(edge_truss_map, truss_order_map, e2, min_e)) {
                    min_e = e2;
                }

                if (test_order(edge_truss_map, truss_order_map, e3, min_e)) {
                    min_e = e3;
                }

                if (min_e == e2) {
                    --rem->at(e2);
                }

                if (min_e == e3) {
                    --rem->at(e3);
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= edge_truss_map->at(e1)) {
                    if (edge_truss_map->at(e2) == min_k) {
                        --ts->at(e2);
                        if (ts->at(e2) < edge_truss_map->at(e2) - 2 && !Q->count(e2)) {
                            Q->insert(e2);
                        }
                    }

                    if (edge_truss_map->at(e3) == min_k) {
                        --ts->at(e3);
                        if (ts->at(e3) < edge_truss_map->at(e3) - 2) {
                            Q->insert(e3);
                        }
                    }
                }
            }

            truss_order_map->at(k)->remove(e1);
            if (!truss_order_map->count(k_new)) {
                truss_order_map->insert(
                        {k_new, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
            }
            truss_order_map->at(k_new)->push_back(e1);
            edge_truss_map->at(e1) = k_new;

            /**
             * @brief recompute ts
             */
            ts->at(e1) = 0;
            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || edge_truss_map->at(e2) < k_new) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_truss_map->at(e3) < k_new) {
                    continue;
                }

                ++ts->at(e1);
            }
        }

//        for (auto p = truss_order_map->at(3)->get_head(); p; p = p->get_next()) {
//            auto e = p->get_value();
//            printf("(%u, %u) ", e->get_source_vertex_id(), e->get_destination_vertex_id());
//        }
//        printf("\n");
    }

    void order_truss_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& current_edge_set,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& C_k,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &s,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                                          const shared_ptr<map<int,shared_ptr<extend_node<int,shared_ptr<abstract_edge>>>>>& B,
                                                          const shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>& e_pivot_node,
                                                          const shared_ptr<abstract_edge>& e_current,
                                                          uint32_t k) {
        auto current_order_list = truss_order_map->at(k);
        auto e_vector = make_shared<vector<shared_ptr<abstract_edge>>>();

        while(!current_edge_set->empty())
        {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set){
                C_k->erase(e1);
                s->erase(e1);

                rem->at(e1) = 0;

                e_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || !(C_k->count(e2) || edge_truss_map->at(e2) > k || current_order_list->count_value(e2) &&
                                                                                    test_order(edge_truss_map,
                                                                                               truss_order_map, e_current,
                                                                                               e2))) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || !(C_k->count(e3) || edge_truss_map->at(e3) > k || current_order_list->count_value(e3) &&
                                                                                 test_order(edge_truss_map, truss_order_map,
                                                                                            e_current, e3))) {
                        continue;
                    }

                    ++rem->at(e1);

                    if (C_k->count(e2)) {
                        --s->at(e2);
                        if (s->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
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
                        if (s->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                            next_edge_set->insert(e3);
                        }
                    } else if (ext->at(e3) > 0) {
                        --ext->at(e3);
                        if (ext->at(e3) == 0) {
                            auto e3_node = current_order_list->find(e3);
                            B->erase(e3_node->get_key());
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*e_vector){
            auto e_node = make_shared<extend_node<int,shared_ptr<abstract_edge>>>(0, e);
            if(!e_pivot_node){
                current_order_list->right_insert(e_node);
            }else{
                current_order_list->insert_before(e_node, e_pivot_node);
            }
        }
    }

    bool order_truss_maintenance::test_order(
            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
            const shared_ptr<unordered_map<uint32_t ,shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
            const shared_ptr<abstract_edge> &e1,
            const shared_ptr<abstract_edge> &e2) {
        auto k1 = edge_truss_map->at(e1);
        auto k2 = edge_truss_map->at(e2);
        if (k1 < k2) {
            return  true;
        }
        else if(k1 > k2){
            return false;
        }else
        {
            auto order_list = truss_order_map->at(k1);
            return order_list->find_key(e1).value() < order_list->find_key(e2).value();
        }
    }

    bool order_truss_maintenance::test_order(const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &order_list,
            const shared_ptr<abstract_edge> &e1,
            const shared_ptr<abstract_edge> &e2) {
            return order_list->find_key(e1).value() < order_list->find_key(e2).value();
    }


    /**
     * @details update the support affected by e
     * @param G
     * @param C_k
     * @param edge_truss_map
     * @param truss_order_map
     * @param ext
     * @param B
     * @param e
     * @param k
     */
    void order_truss_maintenance::update_edge_support(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& C_k,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& ext,
                                                     const shared_ptr<map<int,shared_ptr<extend_node<int,shared_ptr<abstract_edge>>>>>& B,
                                                     const shared_ptr<abstract_edge> &e,
                                                     uint32_t k) {
        auto current_order_list = truss_order_map->at(k);

        const auto &e1 = e;
        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
            swap(u,v);
        }


        for(const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()) {
            if (w == v || !(edge_truss_map->at(e2) > k || current_order_list->count_value(e2)  && test_order(edge_truss_map, truss_order_map, e1, e2))) {
                continue;
            }

            auto e3 = G->get_edge(v, w);
            if (!e3 || !(edge_truss_map->at(e3) > k || current_order_list->count_value(e3)  && test_order(edge_truss_map, truss_order_map, e1, e3))) {
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
        }
    }

    /**
     * @details update the support affected by edges in C_k
     * @param G
     * @param C_k
     * @param edge_truss_map
     * @param truss_order_map
     * @param ext
     * @param B
     * @param e
     * @param k
     */
    void order_truss_maintenance::update_edge_support(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& ext,
                                                     const shared_ptr<map<int,shared_ptr<extend_node<int,shared_ptr<abstract_edge>>>>>& B,
                                                     uint32_t k) {
        auto current_order_list = truss_order_map->at(k);

        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        for(const auto&e1:*C_k){
            visited_edge_set->insert(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                swap(u,v);
            }

            for (const auto &[w, e2]:*G->get_vertex(u)->get_edge_map()) {
                if (w == v || visited_edge_set->count(e2) || !(C_k->count(e2) || edge_truss_map->at(e2) >= k)) {
                    continue;
                }


                auto e3 = G->get_edge(v, w);
                if (!e3 || visited_edge_set->count(e3) || !(C_k->count(e3) || edge_truss_map->at(e3) >= k)) {
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
            }
        }
    }

    void order_truss_maintenance::decrease_rem(const shared_ptr<abstract_graph> &G,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                               const shared_ptr<abstract_edge> &e,
                                               uint32_t k_new) {


        const auto &e1 = e;
        rem->at(e1) = 0;

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

//
//        if(!G->get_vertex(u) || !G->get_vertex(v)){
//            return;
//        }

        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
            swap(u, v);
        }

        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
            if (w == v || edge_truss_map->at(e2) <= k_new) {
                continue;
            }
            auto e3 = G->get_edge(v, w);
            if (!e3 || edge_truss_map->at(e3) <= k_new) {
                continue;
            }

            ++rem->at(e1);

            if (test_order(edge_truss_map, truss_order_map, e2, e1)) {
                --rem->at(e2);
            }

            if (test_order(edge_truss_map, truss_order_map, e3, e1)) {
                --rem->at(e3);
            }
        }
    }

    void order_truss_maintenance::increase_rem(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                              const shared_ptr<abstract_edge> &e) {


        const auto &e1 = e;

        rem->at(e1) = 0;
        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if(!G->get_vertex(u) || !G->get_vertex(v)){
            return;
        }

        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
            swap(u,v);
        }

        for (const auto &[w, e2]:*G->get_vertex(u)->get_edge_map()) {
            if (w == v) {
                continue;
            }
            auto e3 = G->get_edge(v, w);
            if (!e3) {
                continue;
            }

            auto min_e = e1;
            if (test_order(edge_truss_map, truss_order_map, e2, min_e)) {
                min_e = e2;
            }


            if (test_order(edge_truss_map, truss_order_map, e3, min_e)) {
                min_e = e3;
            }

            if (e1 == min_e) {
                ++rem->at(e1);
            }

            if (e2 == min_e) {
                ++rem->at(e2);
            }

            if (e3 == min_e) {
                ++rem->at(e3);
            }
        }
    }

    void order_truss_maintenance::update_rem_and_ts(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> & edge_set,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ts,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Q) {
        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e1:*edge_set) {
            visited_edge_set->insert(e1);

            ts->erase(e1);
            rem->erase(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || visited_edge_set->count(e2)) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || visited_edge_set->count(e3)) {
                    continue;
                }

                if (edge_truss_map->at(e2) > edge_truss_map->at(e1) &&
                    edge_truss_map->at(e3) > edge_truss_map->at(e1)) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (edge_truss_map->at(e2) == min_k) {
                    --ts->at(e2);
                    if (!edge_set->count(e2) && ts->at(e2) < edge_truss_map->at(e2) - 2) {
                        Q->insert(e2);
                    }
                }

                if (edge_truss_map->at(e3) == min_k) {
                    --ts->at(e3);
                    if (!edge_set->count(e3) && ts->at(e3) < edge_truss_map->at(e3) - 2) {
                        Q->insert(e3);
                    }
                }

                auto min_e = e1;

                if (test_order(edge_truss_map, truss_order_map, e2, min_e)) {
                    min_e = e2;
                }

                if (test_order(edge_truss_map, truss_order_map, e3, min_e)) {
                    min_e = e3;
                }

                if (min_e == e2) {
                    --rem->at(e2);
                }

                if (min_e == e3) {
                    --rem->at(e3);
                }
            }
        }
    }


    void order_truss_maintenance::update_order_list(const shared_ptr<abstract_graph> &G,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &s,
                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                   const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                                   const shared_ptr<vector<shared_ptr<abstract_edge>>> &e_vector,
                                                   uint32_t k){
        for(const auto& e1:*C_k)
        {
            s->insert({e1,0});
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                swap(u,v);
            }

            for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v || !C_k->count(e2) && edge_truss_map->at(e2) < k){
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || !C_k->count(e3) && edge_truss_map->at(e3) < k) {
                    continue;
                }

                ++s->at(e1);
            }

            if (s->at(e1) <= k - 2) {
                current_edge_set->insert(e1);
            }
        }

        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set)
            {
                C_k->erase(e1);
                rem->at(e1) = 0;

                e_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u,v);
                }

                for (const auto &[w, e2]:*G->get_vertex(u)->get_edge_map()) {
                    if (w == v || !(C_k->count(e2) || edge_truss_map->at(e2) >= k)) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 ||!(C_k->count(e3) || edge_truss_map->at(e3) >= k)) {
                        continue;
                    }

                    ++rem->at(e1);

                    if (C_k->count(e2)) {
                        --s->at(e2);
                        if (s->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                            next_edge_set->insert(e2);
                        }
                    }

                    if (C_k->count(e3)) {
                        --s->at(e3);
                        if (s->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                            next_edge_set->insert(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }


    void order_truss_maintenance::increase_ts(const shared_ptr<abstract_graph> &G,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &C_k,
                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ts,
                                              uint32_t k) {
        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &e1: *edge_set) {
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            visited_edge_set->insert(e1);

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree())
            {
                swap(u,v);
            }

            for(const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v || visited_edge_set->count(e2) ||!(C_k->count(e2) || edge_truss_map->at(e2) >= k)){
                    continue;
                }
                auto e3 = G->get_edge(v,w);
                if(!e3 || visited_edge_set->count(e3) ||!(C_k->count(e3) || edge_truss_map->at(e3) >= k)){
                    continue;
                }

                ++ts->at(e1);

               if(edge_truss_map->at(e2) == k){
                   ++ts->at(e2);
               }

                if(edge_truss_map->at(e3) == k){
                    ++ts->at(e3);
                }
            }
        }
    }
}