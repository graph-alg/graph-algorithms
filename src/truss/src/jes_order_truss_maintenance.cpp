
#include "truss/jes_order_truss_maintenance.h"

namespace scnu {
    uint32_t jes_order_truss_maintenance::compute_truss_number(const shared_ptr<abstract_graph> &G,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                               const shared_ptr<abstract_edge>& e) {
        auto result_map = make_shared<map<uint32_t,uint32_t> >();
        auto &e1 = e;
        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
            swap(u, v);
        }
        for (const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()) {
            if(w == v){
                continue;
            }
            auto e3 = G->get_edge(v,w);
            if(!e3){
                continue;
            }

            auto min_k = std::min(edge_truss_map->at(e2), edge_truss_map->at(e3));
            if (!result_map->count(min_k)) {
                result_map->insert({min_k, 0});
            }
            ++result_map->at(min_k);
        }

        uint32_t k = 2;
        if(!result_map->empty()){
            uint32_t sum = 0;
            k = result_map->rbegin()->first;
            while (k >= 2) {
                if(result_map->count(k)){
                    sum += result_map->at(k);
                }
                if (sum >= k - 2) {
                    break;
                }
                --k;
            }
        }
        return k;
    }

    void jes_order_truss_maintenance::compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool> &pool) {
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto affected_edge_map = make_shared<unordered_map<shared_ptr<abstract_edge>, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

        auto superior_edge_triangle_map = make_shared<unordered_map<shared_ptr<abstract_edge>, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

        auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(auto iter = ED->begin(); iter!=ED->end(); ){
            auto &e = *iter;
            ++iter;

            if(edge_truss_map->at(e) == 2){
                edge_truss_map->erase(e);
                truss_order_map->at(2)->remove(e);
                G->remove_edge(e);
                ED->erase(e);
            }
            else if (edge_truss_map->at(e) > 2) {
                superior_edge_set->insert(e);
            }
        }

        compute_removal_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);
        /**
         * @brief parallel part
         */
        for (const auto&e1:*Ec) {
            /**
             * @note: thread security
             */
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= edge_truss_map->at(e1)) {
                    if (edge_truss_map->at(e2) == min_k && !Ec->count(e2)) {
                        if (!EU->count(edge_truss_map->at(e2))) {
                            EU->insert(
                                    {edge_truss_map->at(e2), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e2))->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k && !Ec->count(e3)) {
                        if (!EU->count(edge_truss_map->at(e3))) {
                            EU->insert(
                                    {edge_truss_map->at(e3), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e3))->insert(e3);
                    }
                }
            }

            G->remove_edge(e1);
            truss_order_map->at(edge_truss_map->at(e1))->remove(e1);
            ED->erase(e1);
            edge_truss_map->erase(e1);
        }
    }

    void jes_order_truss_maintenance::compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool> &pool) {
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto superior_edge_triangle_map = make_shared<unordered_map<shared_ptr<abstract_edge>, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

        auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(auto iter = ED->begin(); iter!=ED->end(); ){
            auto &e = *iter;
            ++iter;

            if(edge_truss_map->at(e) == 2){
                edge_truss_map->erase(e);
                truss_order_map->at(2)->remove(e);
                G->remove_edge(e);
                ED->erase(e);
            } else if (edge_truss_map->at(e) > 2) {
                superior_edge_set->insert(e);
            }
        }

        compute_removal_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);

        G->insert_edge_collection(Ec);

        for (const auto &e1: *Ec) {
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= edge_truss_map->at(e1)) {
                    if (edge_truss_map->at(e2) == min_k && !Ec->count(e2)) {
                        --TS->at(e2);

                        if (!EU->count(edge_truss_map->at(e2))) {
                            EU->insert(
                                    {edge_truss_map->at(e2), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e2))->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k && !Ec->count(e3)) {
                        --TS->at(e3);

                        if (!EU->count(edge_truss_map->at(e3))) {
                            EU->insert(
                                    {edge_truss_map->at(e3), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e3))->insert(e3);
                    }
                }
            }

            G->remove_edge(e1);
            if (truss_order_map->at(edge_truss_map->at(e1))->count_value(e1)) {
                truss_order_map->at(edge_truss_map->at(e1))->remove(e1);
            }
            ED->erase(e1);
            edge_truss_map->erase(e1);
            TS->erase(e1);
        }
    }

    void jes_order_truss_maintenance::compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool> &pool) {
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(auto iter = ED->begin(); iter!=ED->end(); ){
            auto &e = *iter;
            ++iter;

            if(edge_truss_map->at(e) == 2){
                edge_truss_map->erase(e);
                truss_order_map->at(2)->remove(e);
                G->remove_edge(e);
                ED->erase(e);
            } else if (edge_truss_map->at(e) > 2) {
                superior_edge_set->insert(e);
            }
        }

        compute_removal_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);

        G->insert_edge_collection(Ec);

        for (const auto &e1: *Ec) {
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= edge_truss_map->at(e1)) {
                    if (edge_truss_map->at(e2) == min_k && !Ec->count(e2)) {
                        --TS->at(e2);

                        if (!EU->count(edge_truss_map->at(e2))) {
                            EU->insert(
                                    {edge_truss_map->at(e2), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e2))->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k && !Ec->count(e3)) {
                        --TS->at(e3);

                        if (!EU->count(edge_truss_map->at(e3))) {
                            EU->insert(
                                    {edge_truss_map->at(e3), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e3))->insert(e3);
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

            G->remove_edge(e1);
            truss_order_map->at(edge_truss_map->at(e1))->remove(e1);
            ED->erase(e1);
            edge_truss_map->erase(e1);
            rem->erase(e1);
            TS->erase(e1);
        }
    }


    void jes_order_truss_maintenance::compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>& VE_map,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool> &pool) {
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(auto iter = ED->begin(); iter!=ED->end(); ){
            auto &e = *iter;
            ++iter;

            if(edge_truss_map->at(e) == 2){
                removed_edge_set->insert(e);
                truss_order_map->at(2)->remove(e);
                G->remove_edge(e);
                ED->erase(e);
            }
            else if (edge_truss_map->at(e) > 2) {
                superior_edge_set->insert(e);
            }
        }

        update_VE_map_for_removal(removed_edge_set, edge_truss_map, VE_map);
        for(const auto &e:*removed_edge_set){
            edge_truss_map->erase(e);
        }
        compute_removal_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);
        update_VE_map_for_removal(Ec, edge_truss_map, VE_map);

        for (const auto&e1:*Ec) {
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= edge_truss_map->at(e1)) {
                    if (edge_truss_map->at(e2) == min_k && !Ec->count(e2)) {
                        --TS->at(e2);

                        if (!EU->count(edge_truss_map->at(e2))) {
                            EU->insert(
                                    {edge_truss_map->at(e2), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e2))->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k && !Ec->count(e3)) {
                        --TS->at(e3);

                        if (!EU->count(edge_truss_map->at(e3))) {
                            EU->insert(
                                    {edge_truss_map->at(e3), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                        }
                        EU->at(edge_truss_map->at(e3))->insert(e3);
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

            G->remove_edge(e1);
            truss_order_map->at(edge_truss_map->at(e1))->remove(e1);
            ED->erase(e1);
            edge_truss_map->erase(e1);
            rem->erase(e1);
            TS->erase(e1);
        }
    }

    void jes_order_truss_maintenance::compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool>& pool){
        /**
        * @brief get G+EI
        * @remarks  next loop will operate on current graph
      */
        G->insert_edge_collection(EI);
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        compute_insert_superior_edge_set(G, EI, edge_truss_map, truss_order_map, superior_edge_set, pool);
        compute_insertion_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);

        for(auto iter = EI->begin(); iter!=EI->end();){
            auto &e = *iter;
            ++iter;

            for (const auto &e: *Ec) {
                EI->erase(e);
            }
        }

        /**
         * @brief parallel classify associated edges in Ec
         * @remarks revise edge_truss_map number of vertices in Ec do not affect its neighbors
         */
        for (const auto &e1:*Ec) {
            auto pre_truss_number = compute_pre_truss_number(G, edge_truss_map, e1);
            edge_truss_map->at(e1) = pre_truss_number;

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }
            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if(min_k <= pre_truss_number){
                    if(edge_truss_map->at(e2) == min_k){
                        e1_set->insert(e2);
                    }

                    if(edge_truss_map->at(e3) == min_k){
                        e1_set->insert(e3);
                    }
                }
            }

            for(const auto &e:*e1_set){
                if (!EU->count(edge_truss_map->at(e))) {
                    EU->insert({edge_truss_map->at(e), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                }
                EU->at(edge_truss_map->at(e))->insert(e);

                truss_order_map->at(edge_truss_map->at(e))->remove(e);
            }

            if (!EU->count(pre_truss_number)) {
                EU->insert({pre_truss_number, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }
            EU->at(pre_truss_number)->insert(e1);
        }
    }

    void jes_order_truss_maintenance::compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool>& pool){
        /**
        * @brief get G+EI
        * @remarks  next loop will operate on current graph
      */
        G->insert_edge_collection(EI);
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        compute_insert_superior_edge_set(G, EI, edge_truss_map, truss_order_map, TS, superior_edge_set, EU, pool);
        compute_insertion_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);

        for (const auto &e: *Ec) {
            EI->erase(e);
        }

        /**
         * @brief parallel classify associated edges in Ec
         * @remarks revise edge_truss_map number of vertices in Ec do not affect its neighbors
         */
        for (const auto &e1: *Ec) {
            auto pre_truss_number = compute_pre_truss_number(G, edge_truss_map, e1);
            edge_truss_map->at(e1) = pre_truss_number;

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }
            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if(min_k <= pre_truss_number){
                    if (min_k == pre_truss_number) {
                        ++TS->at(e1);
                    }


                    if (edge_truss_map->at(e2) == min_k) {
                        ++TS->at(e2);
                        e1_set->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k) {
                        ++TS->at(e3);
                        e1_set->insert(e3);
                    }
                } else {
                    ++TS->at(e1);
                }
            }

            for(const auto &e:*e1_set){
                if (!EU->count(edge_truss_map->at(e))) {
                    EU->insert({edge_truss_map->at(e), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                }
                EU->at(edge_truss_map->at(e))->insert(e);

                if(truss_order_map->at(edge_truss_map->at(e))->count_value(e)){
                    truss_order_map->at(edge_truss_map->at(e))->remove(e);
                }
            }

            if (!EU->count(pre_truss_number)) {
                EU->insert({pre_truss_number, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }
            EU->at(pre_truss_number)->insert(e1);
        }
    }

    void jes_order_truss_maintenance::compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool>& pool) {
        /**
        * @brief get G+EI
        * @remarks  next loop will operate on current graph
      */
        G->insert_edge_collection(EI);
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        compute_insert_superior_edge_set(G, EI, edge_truss_map, truss_order_map, rem, TS, superior_edge_set, pool);
        compute_insertion_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);

        for (const auto &e: *Ec) {
            EI->erase(e);
        }

        /**
         * @brief parallel classify associated edges in Ec
         * @remarks revise edge_truss_map number of vertices in Ec do not affect its neighbors
         */
        for (const auto &e1: *Ec) {
            auto pre_truss_number = compute_pre_truss_number(G, edge_truss_map, e1);
            edge_truss_map->at(e1) = pre_truss_number;

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }
            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= pre_truss_number) {
                    if (min_k == pre_truss_number) {
                        ++TS->at(e1);
                    }

                    if (edge_truss_map->at(e2) == min_k) {
                        ++TS->at(e2);
                        e1_set->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k) {
                        ++TS->at(e3);
                        e1_set->insert(e3);
                    }
                } else {
                    ++TS->at(e1);
                }
            }

            for (const auto &e: *e1_set) {
                if (!EU->count(edge_truss_map->at(e))) {
                    EU->insert({edge_truss_map->at(e), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                }
                EU->at(edge_truss_map->at(e))->insert(e);
            }
        }

        for (const auto &[k, k_set]: *EU) {
            pool->submit_task([=] {
                auto current_order_list = truss_order_map->at(k);
                auto visited_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                for (const auto &e1: *k_set) {
                    visited_set->insert(e1);

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();

                    auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                        swap(u, v);
                    }
                    for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                        if (w == v || edge_truss_map->at(e2) < k || visited_set->count(e2) || Ec->count(e2)) {
                            continue;
                        }
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k || visited_set->count(e3) || Ec->count(e3)) {
                            continue;
                        }

                        auto min_e = e1;

                        if (current_order_list->count_value(e2) && test_order(current_order_list, e2, min_e)) {
                            min_e = e2;
                        }

                        if (current_order_list->count_value(e3) && test_order(current_order_list, e3, min_e)) {
                            min_e = e3;
                        }

                        if (min_e == e1) {
                            continue;
                        }

                        if (min_e == e2) {
                            --rem->at(e2);
                        }

                        if (min_e == e3) {
                            --rem->at(e3);
                        }
                    }
                }

                for (const auto &e: *k_set) {
                    truss_order_map->at(k)->remove(e);
                }
            });
        }
        pool->barrier();

        for (const auto &e: *Ec) {
            if (!EU->count(edge_truss_map->at(e))) {
                EU->insert({edge_truss_map->at(e), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }
            EU->at(edge_truss_map->at(e))->insert(e);
        }
    }


    void jes_order_truss_maintenance::compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                              const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                                              const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                              const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU,
                                                              const shared_ptr<thread_pool> &pool) {
        /**
        * @brief get G+EI
        * @remarks  next loop will operate on current graph
      */
        G->insert_edge_collection(EI);
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        compute_insert_superior_edge_set(G, EI, edge_truss_map, truss_order_map, rem, TS, VE_map, superior_edge_set,
                                         pool);
        compute_insertion_maximal_three_hop_independent_set(G, edge_truss_map, superior_edge_set, Ec);

        for (const auto &e: *Ec) {
            EI->erase(e);
        }

        /**
         * @brief parallel classify associated edges in Ec
         * @remarks revise edge_truss_map number of vertices in Ec do not affect its neighbors
         */
        for (const auto &e1: *Ec) {
            auto pre_truss_number = compute_pre_truss_number(G, edge_truss_map, e1);
            edge_truss_map->at(e1) = pre_truss_number;

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }
            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= pre_truss_number) {
                    if (min_k == pre_truss_number) {
                        ++TS->at(e1);
                    }

                    if (edge_truss_map->at(e2) == min_k) {
                        ++TS->at(e2);
                        e1_set->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k) {
                        ++TS->at(e3);
                        e1_set->insert(e3);
                    }
                } else {
                    ++TS->at(e1);
                }
            }

            for (const auto &e: *e1_set) {
                if (!EU->count(edge_truss_map->at(e))) {
                    EU->insert({edge_truss_map->at(e), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                }
                EU->at(edge_truss_map->at(e))->insert(e);
            }
        }

        update_VE_map_for_insertion(Ec,edge_truss_map ,VE_map);

        for (const auto &[k, k_set]: *EU) {
            pool->submit_task([=] {
                auto current_order_list = truss_order_map->at(k);
                auto visited_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
                for (const auto &e1: *k_set) {
                    visited_set->insert(e1);

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();

                    auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                    if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                        swap(u, v);
                    }
                    for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                        if (w == v || edge_truss_map->at(e2) < k || visited_set->count(e2) || Ec->count(e2)) {
                            continue;
                        }
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k || visited_set->count(e3) || Ec->count(e3)) {
                            continue;
                        }

                        auto min_e = e1;

                        if (current_order_list->count_value(e2) && test_order(current_order_list, e2, min_e)) {
                            min_e = e2;
                        }

                        if (current_order_list->count_value(e3) && test_order(current_order_list, e3, min_e)) {
                            min_e = e3;
                        }

                        if (min_e == e1) {
                            continue;
                        }

                        if (min_e == e2) {
                            --rem->at(e2);
                        }

                        if (min_e == e3) {
                            --rem->at(e3);
                        }
                    }
                }

                for (const auto &e: *k_set) {
                    truss_order_map->at(k)->remove(e);
                }
            });
        }
        pool->barrier();

        for (const auto &e: *Ec) {
            if (!EU->count(edge_truss_map->at(e))) {
                EU->insert({edge_truss_map->at(e), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }
            EU->at(edge_truss_map->at(e))->insert(e);
        }
    }

    uint32_t
    jes_order_truss_maintenance::compute_pre_truss_number(const shared_ptr<abstract_graph> &G,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                          const shared_ptr<abstract_edge> &e) {
        return compute_truss_number(G, edge_truss_map, e);
    }


    uint32_t jes_order_truss_maintenance::get_rem_set(const shared_ptr<abstract_graph> &G,
                                                      const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                      const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                      const shared_ptr<abstract_edge>& e,
                                                      const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> & rem_set,
                                                      uint32_t k) {
        auto order_list = truss_order_map->at(k);

        uint32_t rem = 0;
        auto &e1 = e;

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
            swap(u, v);
        }

        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
            if (w == v || edge_truss_map->at(e2) < k || !(edge_truss_map->at(e2) > k || order_list->count_value(e2) &&
                                                                                        test_order(edge_truss_map,
                                                                                                   truss_order_map, e1,
                                                                                                   e2))) {
                continue;
            }


            auto e3 = G->get_edge(v, w);
            if (!e3 || edge_truss_map->at(e3) < k || !(edge_truss_map->at(e3) > k || order_list->count_value(e3) &&
                                                                                     test_order(edge_truss_map,
                                                                                                truss_order_map, e1,
                                                                                                e3))) {
                continue;
            }

            ++rem;

            if(order_list->count_value(e2)){
                rem_set->insert(e2);
            }

            if(order_list->count_value(e3)){
                rem_set->insert(e3);
            }
        }
        return rem;
    }

    shared_ptr<unordered_set<shared_ptr<abstract_edge>>>
    jes_order_truss_maintenance::get_rem_set(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &order_list,
                                             const shared_ptr<abstract_edge> &e,
                                             uint32_t k) {
        auto rem_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        auto &e1 = e;

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
            swap(u, v);
        }

        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
            if (w == v || edge_truss_map->at(e2) < k || !(edge_truss_map->at(e2) > k || order_list->count_value(e2) &&
                                                                                        test_order(order_list, e1,
                                                                                                   e2))) {
                continue;
            }


            auto e3 = G->get_edge(v, w);
            if (!e3 || edge_truss_map->at(e3) < k || !(edge_truss_map->at(e3) > k || order_list->count_value(e3) &&
                                                                                     test_order(order_list, e1,
                                                                                                e3))) {
                continue;
            }

            if(order_list->count_value(e2)){
                rem_set->insert(e2);
            }

            if(order_list->count_value(e3)){
                rem_set->insert(e3);
            }
        }
        return  rem_set;
    }


    shared_ptr<unordered_set<shared_ptr<abstract_edge>>>
    jes_order_truss_maintenance::get_rem_set(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &order_list,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                             const shared_ptr<abstract_edge> &e,
                                             uint32_t k) {
        auto rem_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        auto &e1 = e;

        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
            swap(u, v);
        }

        auto u_map = VE_map->at(u);
        if(u_map->count(k)){
            for(const auto &w:*u_map->at(k)){
                auto e2 = G->get_edge(u, w);
                if (v == w || !order_list->count_value(e2) || test_order(order_list, e2, e1)){
                    continue;
                }

                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_truss_map->at(e3) < k || !(edge_truss_map->at(e3) > k || order_list->count_value(e3) &&
                                                                                         test_order(order_list, e1, e3))) {
                    continue;
                }

                if(order_list->count_value(e2)){
                    rem_set->insert(e2);
                }

                if(order_list->count_value(e3)){
                    rem_set->insert(e3);
                }
            }
        }
        for(auto iter = u_map->lower_bound(k + 1); iter!=u_map->end();++iter){
            for(const auto &w:*iter->second){
                auto e2 = G->get_edge(u, w);
                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_truss_map->at(e3) < k || !(edge_truss_map->at(e3) > k || order_list->count_value(e3) &&
                                                                                         test_order(order_list, e1, e3))) {
                    continue;
                }

                if(order_list->count_value(e3)){
                    rem_set->insert(e3);
                }
            }
        }
        return  rem_set;
    }

    void jes_order_truss_maintenance::compute_insertion_maximal_three_hop_independent_set(
            const shared_ptr<abstract_graph> &G,
            const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
            const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec) {
        G->remove_edge_collection(superior_edge_set);

        auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &e1: *superior_edge_set) {
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            G->insert_edge(e1);

            auto pre_truss_number = compute_pre_truss_number(G, edge_truss_map, e1);

            auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            e1_set->insert(e1);

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v) {
                    continue;
                }

                auto e3 = G->get_edge(v, w);
                if (!e3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= pre_truss_number) {
                    if (edge_truss_map->at(e2) == min_k) {
                        e1_set->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k) {
                        e1_set->insert(e3);
                    }
                }
            }
            bool flag = true;
            for (const auto &e: *e1_set) {
                if (evicted_set->count(e)) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                Ec->insert(e1);
                evicted_set->merge(*e1_set);
            } else {
                G->remove_edge(e1);
            }
        }

    }

    void
    jes_order_truss_maintenance::compute_removal_maximal_three_hop_independent_set(const shared_ptr<abstract_graph> &G,
                                                                                   const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec) {
        for (auto iter = superior_edge_set->begin(); iter != superior_edge_set->end();) {
            auto &e = *iter;
            ++iter;
            if (edge_truss_map->at(e) == 3) {
                Ec->insert(e);
                superior_edge_set->erase(e);
            }
        }

        auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &e1: *superior_edge_set) {

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            e1_set->insert(e1);

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || edge_truss_map->at(e2) == 3) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_truss_map->at(e3) == 3) {
                    continue;
                }

                auto min_k = min(edge_truss_map->at(e2), edge_truss_map->at(e3));
                if (min_k <= edge_truss_map->at(e1)) {
                    if (edge_truss_map->at(e2) == min_k) {
                        e1_set->insert(e2);
                    }

                    if (edge_truss_map->at(e3) == min_k) {
                        e1_set->insert(e3);
                    }
                }
            }

            bool flag = true;
            for (const auto &e: *e1_set) {
                if (evicted_set->count(e)) {
                    flag = false;
                }
            }
            if (flag) {
                Ec->insert(e1);
                evicted_set->merge(*e1_set);
                G->remove_edge(e1);
            }
        }
    }

    uint32_t jes_order_truss_maintenance::compute_truss_support(const shared_ptr<abstract_graph> &G,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                const shared_ptr<abstract_edge>& e,
                                                                uint32_t k) {
        uint32_t edge_support = 0;
        auto &e1 = e;
        auto u = e1->get_source_vertex_id();
        auto v = e1->get_destination_vertex_id();

        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree())
        {
            swap(u, v);
        }

        for (const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()) {
            if(w == v || edge_truss_map->at(e2) < k)
            {
                continue;
            }
            auto e3 = G->get_edge(v,w);
            if(!e3 || edge_truss_map->at(e3) < k)
            {
                continue;
            }

            ++edge_support;
        }
        return edge_support;
    }


    void jes_order_truss_maintenance::compute_insert_superior_edge_set(const shared_ptr<abstract_graph>&G,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& superior_edge_set,
                                                                       const shared_ptr<thread_pool>& pool) {
        auto edge_truss_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        for (const auto &e:*edge_set) {
            edge_truss_support_map->insert({e, 0});
        }

        auto thread_number = pool->get_thread_number();
        {
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{

                    for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                            swap(u, v);
                        }

                        for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                            if(w == v){
                                continue;
                            }
                            auto e3 = G->get_edge(v,w);
                            if(!e3 ){
                                continue;
                            }

                            ++edge_truss_support_map->at(e1);
                        }
                    }
                });
            }
        }
        pool->barrier();

        if (!truss_order_map->count(2)) {
            truss_order_map->insert({2, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
        }

        auto sub_inserted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &[e, e_support]: *edge_truss_support_map) {
            if (e_support == 0) {
                edge_set->erase(e);
                truss_order_map->at(2)->push_back(e);
            }
            if (e_support == 1) {
                edge_set->erase(e);
                sub_inserted_edge_set->insert(e);
            }

            if (e_support > 1) {
                superior_edge_set->insert(e);
            }
        }

        if (!sub_inserted_edge_set->empty()) {
            auto TS = make_shared<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>>();
            for (const auto &e: *sub_inserted_edge_set) {
                TS->insert({e, 0});
            }

            /**
             * @brief inserted edges with support 1
             */
            auto candidate_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            auto visited_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            for (const auto &e1: *sub_inserted_edge_set) {
                visited_edge_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || superior_edge_set->count(e2) || visited_edge_set->count(e2)) {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if (!e3 || superior_edge_set->count(e3) || visited_edge_set->count(e3)) {
                        continue;
                    }

                    ++TS->at(e1);

                    if (edge_truss_map->at(e2) == 2) {
                        if (!sub_inserted_edge_set->count(e2)) {
                            candidate_edge_set->insert(e2);
                        } else {
                            ++TS->at(e2);
                        }
                    }

                    if (edge_truss_map->at(e3) == 2) {
                        if (!sub_inserted_edge_set->count(e3)) {
                            candidate_edge_set->insert(e3);
                        } else {
                            ++TS->at(e3);
                        }
                    }
                }
            }

            for (auto iter = sub_inserted_edge_set->begin(); iter != sub_inserted_edge_set->end();) {
                auto &e = *iter;
                ++iter;

                if (edge_truss_support_map->at(e) == 0) {
                    truss_order_map->at(2)->push_back(e);
                    sub_inserted_edge_set->erase(e);
                }
            }

            for (const auto &e: *candidate_edge_set) {
                truss_order_map->at(2)->remove(e);
            }
            candidate_edge_set->merge(*sub_inserted_edge_set);

            for (const auto &e: *candidate_edge_set) {
                edge_truss_map->at(e) = 3;
            }

            if (!candidate_edge_set->empty()) {
                if (!truss_order_map->count(3)) {
                    truss_order_map->insert({3, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                }

                auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                k_joint_insert(G, candidate_edge_set, edge_truss_map, candidate_edge_support_map, truss_order_map, 3);
            }
        }
    }

    void jes_order_truss_maintenance::compute_insert_superior_edge_set(const shared_ptr<abstract_graph>&G,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& superior_edge_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>& EU,
                                                                       const shared_ptr<thread_pool>& pool) {
        auto edge_truss_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        for (const auto &e:*edge_set) {
            edge_truss_support_map->insert({e, 0});
        }

        auto thread_number = pool->get_thread_number();
        {
            auto location_vector = pool->split_task(edge_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{

                    for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                            swap(u, v);
                        }

                        for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                            if(w == v){
                                continue;
                            }
                            auto e3 = G->get_edge(v,w);
                            if(!e3 ){
                                continue;
                            }

                            ++edge_truss_support_map->at(e1);
                        }
                    }
                });
            }
        }
        pool->barrier();

        if (!truss_order_map->count(2)) {
            truss_order_map->insert({2, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
        }

        auto sub_inserted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &[e, e_support]: *edge_truss_support_map) {
            if (e_support == 0) {
                edge_set->erase(e);
                truss_order_map->at(2)->push_back(e);
            }
            if (e_support == 1) {
                edge_set->erase(e);
                sub_inserted_edge_set->insert(e);
            }

            if (e_support > 1) {
                superior_edge_set->insert(e);
            }
        }

        if (!sub_inserted_edge_set->empty()) {
            /**
             * @brief inserted edges with support 1
             */
            auto candidate_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            auto visited_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            for (const auto &e1: *sub_inserted_edge_set) {
                visited_edge_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || superior_edge_set->count(e2) || visited_edge_set->count(e2)) {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if (!e3 || superior_edge_set->count(e3) || visited_edge_set->count(e3)) {
                        continue;
                    }

                    ++TS->at(e1);

                    if (edge_truss_map->at(e2) == 2) {
                        ++TS->at(e2);
                        if (!sub_inserted_edge_set->count(e2)) {
                            candidate_edge_set->insert(e2);
                        }
                    }

                    if (edge_truss_map->at(e3) == 2) {
                        ++TS->at(e3);
                        if (!sub_inserted_edge_set->count(e3)) {
                            candidate_edge_set->insert(e3);
                        }
                    }
                }
            }

            for (auto iter = sub_inserted_edge_set->begin(); iter != sub_inserted_edge_set->end();) {
                auto &e = *iter;
                ++iter;

                if (TS->at(e) == 0) {
                    truss_order_map->at(2)->push_back(e);
                    sub_inserted_edge_set->erase(e);
                }
            }

            for (const auto &e: *candidate_edge_set) {
                truss_order_map->at(2)->remove(e);
            }
            candidate_edge_set->merge(*sub_inserted_edge_set);

            for (const auto &e: *candidate_edge_set) {
                edge_truss_map->at(e) = 3;
            }

            if (!candidate_edge_set->empty()) {
                update_edge_truss_support(G, candidate_edge_set, edge_truss_map, TS, 3);

                if(!truss_order_map->count(3)){
                    truss_order_map->insert({3, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                }
//                auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
//                k_joint_insert(G, candidate_edge_set, edge_truss_map, candidate_edge_support_map,
//                                       truss_order_map, TS, 3);
//
                if(!EU->count(3)){
                    EU->insert({3, candidate_edge_set});
                }else{
                    EU->at(3)->merge(*candidate_edge_set);
                }
            }
        }
    }

    void jes_order_truss_maintenance::compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                       const shared_ptr<thread_pool> &pool) {
        auto edge_truss_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        for (const auto &e: *edge_set) {
            edge_truss_support_map->insert({e, 0});
        }

        auto thread_number = pool->get_thread_number();
        {
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {

                    for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                            swap(u, v);
                        }

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v) {
                                continue;
                            }
                            auto e3 = G->get_edge(v, w);
                            if (!e3) {
                                continue;
                            }

                            ++edge_truss_support_map->at(e1);
                        }
                    }
                });
            }
        }
        pool->barrier();

        if (!truss_order_map->count(2)) {
            truss_order_map->insert({2, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
        }

        for (const auto &[e, e_support]: *edge_truss_support_map) {
            if (e_support == 0) {
                edge_set->erase(e);
                truss_order_map->at(2)->push_back(e);
            }

            if (e_support >= 1) {
                superior_edge_set->insert(e);
            }
        }
    }


    void jes_order_truss_maintenance::compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                       const shared_ptr<thread_pool> &pool) {
        auto edge_truss_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        for (const auto &e: *edge_set) {
            edge_truss_support_map->insert({e, 0});
        }

        auto thread_number = pool->get_thread_number();
        {
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {

                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                            swap(u, v);
                        }

                        for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                            if(w == v){
                                continue;
                            }
                            auto e3 = G->get_edge(v,w);
                            if(!e3 ){
                                continue;
                            }

                            ++edge_truss_support_map->at(e1);
                        }
                    }
                });
            }
        }
        pool->barrier();

        if (!truss_order_map->count(2)) {
            truss_order_map->insert({2, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
        }

        auto sub_inserted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &[e, e_support]: *edge_truss_support_map) {
            if (e_support == 0) {
                edge_set->erase(e);
                truss_order_map->at(2)->push_back(e);
            }
            if (e_support == 1) {
                edge_set->erase(e);
                sub_inserted_edge_set->insert(e);
            }

            if (e_support > 1) {
                superior_edge_set->insert(e);
            }
        }

        if (!sub_inserted_edge_set->empty()) {
            /**
             * @brief inserted edges with support 1
             */
            auto candidate_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            auto visited_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            for (const auto &e1: *sub_inserted_edge_set) {
                visited_edge_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || superior_edge_set->count(e2) || visited_edge_set->count(e2)) {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if (!e3 || superior_edge_set->count(e3) || visited_edge_set->count(e3)) {
                        continue;
                    }

                    ++TS->at(e1);

                    if (edge_truss_map->at(e2) == 2) {
                        ++TS->at(e2);
                        if (!sub_inserted_edge_set->count(e2)) {
                            candidate_edge_set->insert(e2);
                        }
                    }

                    if (edge_truss_map->at(e3) == 2) {
                        ++TS->at(e3);
                        if (!sub_inserted_edge_set->count(e3)) {
                            candidate_edge_set->insert(e3);
                        }
                    }
                }
            }

            for (auto iter = sub_inserted_edge_set->begin(); iter != sub_inserted_edge_set->end();) {
                auto &e = *iter;
                ++iter;

                if (TS->at(e) == 0) {
                    truss_order_map->at(2)->push_back(e);
                    sub_inserted_edge_set->erase(e);
                }
            }

            for (const auto &e: *candidate_edge_set) {
                truss_order_map->at(2)->remove(e);
            }
            candidate_edge_set->merge(*sub_inserted_edge_set);

            for (const auto &e: *candidate_edge_set) {
                edge_truss_map->at(e) = 3;
            }

            if (!candidate_edge_set->empty()) {
                update_edge_truss_support(G, candidate_edge_set, edge_truss_map, TS, 3);

                if (!truss_order_map->count(3)) {
                    truss_order_map->insert({3, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                }

                auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                k_joint_insert(G, candidate_edge_set, edge_truss_map, candidate_edge_support_map, truss_order_map, rem,
                               TS,
                               3);
            }
        }
    }


    void jes_order_truss_maintenance::compute_insert_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                       const shared_ptr<thread_pool> &pool) {
        auto edge_truss_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
        for (const auto &e: *edge_set) {
            edge_truss_support_map->insert({e, 0});
        }

        auto thread_number = pool->get_thread_number();
        {
            auto location_vector = pool->split_task(edge_set);
            for (uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {

                    for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                        auto &e1 = *iter;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                            swap(u, v);
                        }

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v) {
                                continue;
                            }
                            auto e3 = G->get_edge(v, w);
                            if (!e3) {
                                continue;
                            }

                            ++edge_truss_support_map->at(e1);
                        }
                    }
                });
            }
        }
        pool->barrier();

        if (!truss_order_map->count(2)) {
            truss_order_map->insert({2, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
        }

        auto direct_inserted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        auto sub_inserted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &[e, e_support]: *edge_truss_support_map) {
            if (e_support == 0) {
                edge_set->erase(e);
                direct_inserted_edge_set->insert(e);
                truss_order_map->at(2)->push_back(e);
            }
            if (e_support == 1) {
                edge_set->erase(e);
                sub_inserted_edge_set->insert(e);
            }

            if (e_support > 1) {
                superior_edge_set->insert(e);
            }
        }

        if (!sub_inserted_edge_set->empty()) {
            /**
             * @brief inserted edges with support 1
             */
            auto candidate_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            auto visited_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
            for (const auto &e1: *sub_inserted_edge_set) {
                visited_edge_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || superior_edge_set->count(e2) || visited_edge_set->count(e2)) {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if (!e3 || superior_edge_set->count(e3) || visited_edge_set->count(e3)) {
                        continue;
                    }

                    ++TS->at(e1);

                    if (edge_truss_map->at(e2) == 2) {
                        ++TS->at(e2);
                        if (!sub_inserted_edge_set->count(e2)) {
                            candidate_edge_set->insert(e2);
                        }
                    }

                    if (edge_truss_map->at(e3) == 2) {
                        ++TS->at(e3);
                        if (!sub_inserted_edge_set->count(e3)) {
                            candidate_edge_set->insert(e3);
                        }
                    }
                }
            }

            for (auto iter = sub_inserted_edge_set->begin(); iter != sub_inserted_edge_set->end();) {
                auto &e = *iter;
                ++iter;

                if (TS->at(e) == 0) {
                    direct_inserted_edge_set->insert(e);
                    truss_order_map->at(2)->push_back(e);
                    sub_inserted_edge_set->erase(e);
                }
            }

            update_VE_map_for_insertion(direct_inserted_edge_set, edge_truss_map, VE_map);

            for (const auto &e: *candidate_edge_set) {
                truss_order_map->at(2)->remove(e);
            }
            update_VE_map_for_removal(candidate_edge_set, edge_truss_map, VE_map);

            candidate_edge_set->merge(*sub_inserted_edge_set);
            for (const auto &e: *candidate_edge_set) {
                edge_truss_map->at(e) = 3;
            }

            if (!candidate_edge_set->empty()) {
                update_VE_map_for_insertion(candidate_edge_set, edge_truss_map, VE_map);

                update_edge_truss_support(G, candidate_edge_set, edge_truss_map, TS, VE_map, 3);

                if (!truss_order_map->count(3)) {
                    truss_order_map->insert({3, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                }

                auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                k_joint_insert(G, candidate_edge_set, edge_truss_map, candidate_edge_support_map, truss_order_map, rem,
                               TS,
                               3);
            }
        }
    }

    void jes_order_truss_maintenance::compute_delete_superior_edge_set(const shared_ptr<abstract_graph> &G,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &affected_edge_map,
                                                                       const shared_ptr<thread_pool> &pool) {
        auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (auto iter = edge_set->begin(); iter != edge_set->end();) {
            auto &e = *iter;
            ++iter;

            if(edge_truss_map->at(e) == 2){
                edge_truss_map->erase(e);
                G->remove_edge(e);
                edge_set->erase(e);
            }
//            else if (edge_truss_map->at(e) == 3) {
//                superior_edge_set->insert(e);
//            }
            else if (edge_truss_map->at(e) > 2) {
                superior_edge_set->insert(e);
            }
        }

//        auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
//        for (const auto &e1: *removed_edge_set) {
//
//            auto u = e1->get_source_vertex_id();
//            auto v = e1->get_destination_vertex_id();
//
//            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
//                swap(u, v);
//            }
//
//            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
//                if (w == v) {
//                    continue;
//                }
//                auto e3 = G->get_edge(v, w);
//                if (!e3) {
//                    continue;
//                }
//
//                if (edge_truss_map->at(e2) == 3 && !removed_edge_set->count(e2)) {
//                    candidate_edge_set->insert(e2);
//                }
//
//                if (edge_truss_map->at(e3) == 3 && !removed_edge_set->count(e3)) {
//                    candidate_edge_set->insert(e3);
//                }
//            }
//
//            edge_set->erase(e1);
//            edge_truss_map->erase(e1);
//            G->remove_edge(e1);
//        }
//
//        if (!candidate_edge_set->empty()) {
//            if (!truss_order_map->count(2)) {
//                truss_order_map->insert({2, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
//            }
//
//            auto candidate_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
//
//            k_joint_delete(G, candidate_edge_set, edge_truss_map, candidate_edge_vector, 3);
//
//
//            for (const auto &e: *candidate_edge_vector) {
//                truss_order_map->at(3)->remove(e);
//                edge_truss_map->at(e) = 2;
//                truss_order_map->at(2)->push_back(e);
//            }
//        }
    }

    void jes_order_truss_maintenance::compute_delete_superior_edge_set(const shared_ptr<abstract_graph>&G,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                       const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& superior_edge_set,
                                                                       const shared_ptr<thread_pool>& pool) {
        auto result_map = make_shared<unordered_map<shared_ptr<abstract_edge>,uint32_t>>();
        for (const auto &e:*edge_set) {
            result_map->insert({e, 0});
        }

        auto location_vector = pool->split_task(edge_set);
        for(uint32_t i = 0; i < pool->get_thread_number(); ++i){
            pool->submit_task([=]{
                for(auto iter = *location_vector->at(i); iter != *location_vector->at(i+1); ++iter){
                    auto &e1 = *iter;

                    auto u = e1->get_source_vertex_id();
                    auto v = e1->get_destination_vertex_id();

                    if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                        swap(u, v);
                    }

                    for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                        if(w == v){
                            continue;
                        }
                        auto e3 = G->get_edge(v,w);
                        if(!e3){
                            continue;
                        }

                        ++result_map->at(e1);
                    }
                }
            });
        }
        pool->barrier();

        auto removed_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto& [e,superior_support]:*result_map)
        {
            if(superior_support == 0){
                edge_set->erase(e);
                G->remove_edge(e);
                edge_truss_map->erase(e);
                TS->erase(e);
            }

            if(superior_support == 1){
                removed_edge_set->insert(e);
            }
            if(superior_support >= 2){
                superior_edge_set->insert(e);
            }
        }

        auto candidate_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e1:*removed_edge_set){

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                swap(u, v);
            }

            for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v){
                    continue;
                }
                auto e3 = G->get_edge(v,w);
                if(!e3){
                    continue;
                }

                if(edge_truss_map->at(e2) == 3 && !removed_edge_set->count(e2)){
                    --TS->at(e2);
                    candidate_edge_set->insert(e2);
                }

                if(edge_truss_map->at(e3) == 3 && !removed_edge_set->count(e3)){
                    --TS->at(e3);
                    candidate_edge_set->insert(e3);
                }
            }

            edge_set->erase(e1);
            edge_truss_map->erase(e1);
            TS->erase(e1);
            G->remove_edge(e1);
        }

        if(!candidate_edge_set->empty()){
            if(!truss_order_map->count(2)){
                truss_order_map->insert({2, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
            }

            auto candidate_edge_vector = make_shared<vector<shared_ptr<abstract_edge>>>();
            k_joint_delete(G, candidate_edge_set, edge_truss_map, TS, candidate_edge_vector, 3);
            for(const auto &e:*candidate_edge_vector){
                truss_order_map->at(3)->remove(e);
                edge_truss_map->at(e) = 2;
                TS->at(e) = 0;
                truss_order_map->at(2)->push_back(e);
            }
        }
    }


    void jes_order_truss_maintenance::edge_support_computation(const shared_ptr<abstract_graph> &G,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>,uint32_t>>& candidate_edge_support_map,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &evicted_edge_set,
                                                               uint32_t k){
        for (const auto &e1: *candidate_edge_set) {
            candidate_edge_support_map->insert({e1, 0});
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || edge_truss_map->at(e2) < k) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_truss_map->at(e3) < k) {
                    continue;
                }

                ++candidate_edge_support_map->at(e1);
            }

            if (candidate_edge_support_map->at(e1) <= k - 2) {
                evicted_edge_set->insert(e1);
            }
        }
    }

    void jes_order_truss_maintenance::init(const shared_ptr<abstract_graph> &G,
                                           const std::shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                           const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                           const shared_ptr<thread_pool> &pool) {
        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            VE_map->insert({u, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
        }

        for (const auto &[u, u_vertex]: *G->get_vertex_map()) {
            pool->submit_task([=] {
                VE_map->at(u) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                auto u_map = VE_map->at(u);
                for (const auto &[v, e]: *u_vertex->get_edge_map()) {
                    auto k = edge_truss_map->at(e);
                    if (!u_map->count(k)) {
                        u_map->insert({k, make_shared<unordered_set<uint32_t>>()});
                    }
                    u_map->at(k)->insert(v);
                }
            });
        }
        pool->barrier();
    }

    void jes_order_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<thread_pool> &pool) {

        for (const auto &e: *EI) {
            edge_truss_map->insert({e, 2});
        }

        while (!EI->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_insert_edge_set(G, EI, edge_truss_map, truss_order_map, Ec, EU, pool);

            auto candidate_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            for(const auto &[k, k_set]:*EU){
                candidate_edge_map->insert({k, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>()});
            }

            /**
              * @brief parallel part
              */
            for (const auto &k_pair: *EU) {
                pool->submit_task([=] {
                    auto k = k_pair.first;
                    auto EUk = k_pair.second;
                    candidate_edge_map->at(k) = make_shared<unordered_set<shared_ptr<abstract_edge>>>(*EUk);
                    auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                    k_joint_insert(G, candidate_edge_map->at(k), edge_truss_map, candidate_edge_support_map,
                                   truss_order_map, k);
                });
            }
            pool->barrier();

            for (const auto &[k, e_set]: *candidate_edge_map) {
                for(const auto &e:*e_set){
                    edge_truss_map->at(e) =  k + 1;
                }
            }

            for(const auto&[k, e_set]:*candidate_edge_map) {
                if (!truss_order_map->count(k + 1)) {
                    truss_order_map->insert({k + 1, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                }
                pool->submit_task([=] {
                    auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                    k_joint_insert(G, e_set, edge_truss_map, candidate_edge_support_map,
                                   truss_order_map, k + 1);
                });
            }
           pool->barrier();
        }
    }

    uint32_t jes_order_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                 const shared_ptr<thread_pool> &pool) {

        for (const auto &e: *EI) {
            edge_truss_map->insert({e, 2});
            TS->insert({e, 0});
        }

        auto next_EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

        uint32_t loop_count = 0;
        while (!EI->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_insert_edge_set(G, EI, edge_truss_map, truss_order_map, TS, Ec, EU, pool);

            ++loop_count;

            for(const auto &[k, e_set]:*next_EU){
                 if(!EU->count(k)){
                     EU->insert({k, e_set});
                 }else{
                     EU->at(k)->merge(*e_set);
                 }
            }

            auto candidate_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            auto candidate_support_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>>>();
            for (const auto &[k, k_set]: *EU) {
                candidate_edge_map->insert({k, k_set});
                candidate_support_map->insert({k, shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>()});
            }

            /**
              * @brief parallel part
              */
            for (const auto &[k, EUk]: *EU) {
                pool->submit_task([=] {
                    candidate_support_map->at(k) = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                    k_joint_insert(G, candidate_edge_map->at(k), edge_truss_map, candidate_support_map->at(k),
                                   truss_order_map, TS, k);
                });
            }
            pool->barrier();

            for (const auto &[k, e_set]: *candidate_edge_map) {
                if (!e_set->empty()) {
                    if (!truss_order_map->count(k + 1)) {
                        truss_order_map->insert({k + 1, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                    }
                    pool->submit_task([=] {
                        for (const auto &e: *e_set) {
                            edge_truss_map->at(e) = k + 1;
                            auto edge_truss_support_map = candidate_support_map->at(k);
                            TS->at(e) = edge_truss_support_map->at(e);
                        }
                    });
                }
            }
            next_EU->clear();
            pool->barrier();

            for (const auto &[k, e_set]: *candidate_edge_map) {
                if (!e_set->empty()) {
                    next_EU->insert({k + 1, e_set});
                    pool->submit_task([=] {
                        update_edge_truss_support(G, e_set, edge_truss_map, TS, k + 1);
//                        auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
//                        k_joint_insert(G, e_set, edge_truss_map, candidate_edge_support_map,
//                                       truss_order_map, TS, k + 1);
                    });
                }
            }
            pool->barrier();
        }
        return loop_count;
    }

    uint32_t jes_order_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                 const shared_ptr<thread_pool> &pool) {

        for (const auto &e: *EI) {
            edge_truss_map->insert({e, 2});
            TS->insert({e, 0});
            rem->insert({e, 0});
        }

        uint32_t loop_count = 0;
        while (!EI->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_insert_edge_set(G, EI, edge_truss_map, truss_order_map, rem, TS, Ec, EU, pool);

            ++loop_count;

            auto candidate_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            auto candidate_support_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>>>();
            for (const auto &[k, k_set]: *EU) {
                candidate_edge_map->insert({k, k_set});
                candidate_support_map->insert({k, shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>()});
            }

            /**
              * @brief parallel part
              */
            for (const auto &[k, EUk]: *EU) {
                pool->submit_task([=] {
                    candidate_support_map->at(k) = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                    k_joint_insert(G, candidate_edge_map->at(k), edge_truss_map, candidate_support_map->at(k),
                                   truss_order_map, rem, TS, k);
                });
            }
            pool->barrier();

            for (const auto &[k, e_set]: *candidate_edge_map) {
                if (!e_set->empty()) {
                    if (!truss_order_map->count(k + 1)) {
                        truss_order_map->insert({k + 1, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                    }
                    pool->submit_task([=] {
                        for (const auto &e: *e_set) {
                            edge_truss_map->at(e) = k + 1;
                            auto edge_truss_support_map = candidate_support_map->at(k);
                            TS->at(e) = edge_truss_support_map->at(e);
                        }
                    });
                }
            }
            pool->barrier();

            for (const auto &[k, e_set]: *candidate_edge_map) {
                if (!e_set->empty()) {
                    pool->submit_task([=] {
                        update_edge_truss_support(G, e_set, edge_truss_map, TS, k + 1);
                        auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                        k_joint_insert(G, e_set, edge_truss_map, candidate_edge_support_map,
                                       truss_order_map, rem, TS, k + 1);
                    });
                }
            }
            pool->barrier();
        }

        return loop_count;
    }


    void jes_order_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                             const shared_ptr<thread_pool> &pool) {

        for (const auto &e: *EI) {
            edge_truss_map->insert({e, 2});
            TS->insert({e, 0});
            rem->insert({e, 0});

            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
        }

        while (!EI->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_insert_edge_set(G, EI, edge_truss_map, truss_order_map, rem, TS, VE_map, Ec, EU, pool);

            auto candidate_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            auto candidate_support_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>>>();
            for (const auto &[k, k_set]: *EU) {
                candidate_edge_map->insert({k, k_set});
                candidate_support_map->insert({k, shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>()});
            }

            /**
              * @brief parallel part
              */
            for (const auto &[k, EUk]: *EU) {
                pool->submit_task([=] {
                    candidate_support_map->at(k) = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                    k_joint_insert(G, candidate_edge_map->at(k), edge_truss_map, candidate_support_map->at(k),
                                   truss_order_map, rem, TS, VE_map, k);
                });
            }
            pool->barrier();


            for (auto iter = candidate_edge_map->begin(); iter != candidate_edge_map->end();) {
                auto [k, e_set] = *iter;
                ++iter;
                if (!e_set->empty()) {
                    update_VE_map_for_removal(e_set, edge_truss_map, VE_map);
                } else {
                    candidate_edge_map->erase(k);
                }
            }

            for (const auto &[k, e_set]: *candidate_edge_map) {
                pool->submit_task([=] {
                    for (const auto &e: *e_set) {
                        edge_truss_map->at(e) = k + 1;
                        auto edge_truss_support_map = candidate_support_map->at(k);
                        TS->at(e) = edge_truss_support_map->at(e);
                    }
                });
            }
            pool->barrier();

            for (const auto &[k, e_set]: *candidate_edge_map) {
                if (!e_set->empty()) {
                    update_VE_map_for_insertion(e_set, edge_truss_map, VE_map);
                }
            }

            for (const auto &[k, e_set]: *candidate_edge_map) {
                if (!e_set->empty()) {
                    if (!truss_order_map->count(k + 1)) {
                        truss_order_map->insert({k + 1, make_shared<extend_list<int, shared_ptr<abstract_edge>>>()});
                    }
                    pool->submit_task([=] {
                        update_edge_truss_support(G, e_set, edge_truss_map, TS, VE_map, k + 1);
                        auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                        k_joint_insert(G, e_set, edge_truss_map, candidate_edge_support_map,
                                       truss_order_map, rem, TS, VE_map, k + 1);
                    });
                }
            }
            pool->barrier();
        }
    }


    void jes_order_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                             const shared_ptr<thread_pool> &pool) {
        while (!ED->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_delete_edge_set(G, ED, edge_truss_map, truss_order_map, Ec, EU, pool);

            auto candidate_edge_map = make_shared<map<uint32_t,shared_ptr<vector<shared_ptr<abstract_edge>>>>>();
            for (const auto &[k,EUk]:*EU) {
                candidate_edge_map->insert({k, make_shared<vector<shared_ptr<abstract_edge>>>()});
            }

            /**
             * @brief parallel part
             */
            for (const auto &k_pair:*EU) {
                pool->submit_task([=]() {
                    auto k = k_pair.first;
                    auto Ek = k_pair.second;
                    k_joint_delete(G, Ek, edge_truss_map,candidate_edge_map->at(k), k);
                });
            }
            pool->barrier();

            for(const auto&[k, k_vector]:*candidate_edge_map)
            {
                for (const auto &e: *k_vector) {
                    truss_order_map->at(k)->remove(e);
                    --edge_truss_map->at(e);
                    truss_order_map->at(k - 1)->push_back(e);
                }
            }
            pool->barrier();
        }
    }

    uint32_t jes_order_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                 const shared_ptr<thread_pool> &pool) {
        uint32_t loop_count = 0;
        while (!ED->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_delete_edge_set(G, ED, edge_truss_map, truss_order_map, TS, Ec, EU, pool);

            ++loop_count;

            auto candidate_edge_map = make_shared<map<uint32_t, shared_ptr<vector<shared_ptr<abstract_edge>>>>>();
            for (const auto &[k, EUk]: *EU) {
                candidate_edge_map->insert({k, make_shared<vector<shared_ptr<abstract_edge>>>()});
            }

            /**
             * @brief parallel part
             */
            for (const auto &k_pair: *EU) {
                pool->submit_task([=]() {
                    auto k = k_pair.first;
                    auto Ek = k_pair.second;
                    k_joint_delete(G, Ek, edge_truss_map, TS, candidate_edge_map->at(k), k);
                });
            }
            pool->barrier();


            for(const auto&[k, k_vector]:*candidate_edge_map) {
                for (const auto &e: *k_vector) {
                    if (truss_order_map->at(k)->count_value(e)) {
                        truss_order_map->at(k)->remove(e);
                    }
                    --edge_truss_map->at(e);
                    truss_order_map->at(k - 1)->push_back(e);
                }
            }

            for (const auto &[k, k_vector]: *candidate_edge_map) {
                pool->submit_task([=] {
                    for (const auto &e: *k_vector) {
                        TS->at(e) = compute_truss_support(G, edge_truss_map, e, k - 1);
                    }
                });
            }
            pool->barrier();
        }

        return loop_count;
    }

    uint32_t jes_order_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                 const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                 const shared_ptr<thread_pool> &pool) {
        uint32_t loop_count = 0;
        while (!ED->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_delete_edge_set(G, ED, edge_truss_map, truss_order_map, rem, TS, Ec, EU, pool);

            ++loop_count;

            auto candidate_edge_map = make_shared<map<uint32_t, shared_ptr<vector<shared_ptr<abstract_edge>>>>>();
            for (const auto &[k, EUk]: *EU) {
                candidate_edge_map->insert({k, make_shared<vector<shared_ptr<abstract_edge>>>()});
            }

            /**
             * @brief parallel part
             */
            for (const auto &k_pair: *EU) {
                pool->submit_task([=]() {
                    auto k = k_pair.first;
                    auto Ek = k_pair.second;
                    k_joint_delete(G, Ek, edge_truss_map, truss_order_map, rem, TS, candidate_edge_map->at(k), k);
                });
            }
            pool->barrier();


            for(const auto&[k, k_vector]:*candidate_edge_map) {
                pool->submit_task([=]{
                    for (const auto &e: *k_vector) {
                        truss_order_map->at(k)->remove(e);
                        --edge_truss_map->at(e);
                    }
                });
            }
            pool->barrier();

            for(const auto&[k, k_vector]:*candidate_edge_map) {
                pool->submit_task([=]{
                    for (const auto &e: *k_vector) {
                        truss_order_map->at(k - 1)->push_back(e);
                    }
                });
            }
            pool->barrier();

            for(const auto&[k, k_vector]:*candidate_edge_map)
            {
                pool->submit_task([=] {
                for (const auto &e1: *k_vector) {
                        TS->at(e1) = 0;
                        rem->at(e1) = 0;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                            swap(u, v);
                        }

                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                            if (w == v || edge_truss_map->at(e2) < k - 1) {
                                continue;
                            }
                            auto e3 = G->get_edge(v, w);
                            if (!e3 || edge_truss_map->at(e3) <  k - 1) {
                                continue;
                            }

                            ++TS->at(e1);

                            if (test_order(edge_truss_map, truss_order_map, e1, e2) &&
                                test_order(edge_truss_map, truss_order_map, e1, e3)) {
                                ++rem->at(e1);
                            }
                        }
                    }
                });
            }
            pool->barrier();
        }
        return loop_count;
    }

    void jes_order_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                             const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int , shared_ptr<abstract_edge>>>>>& truss_order_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& rem,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& TS,
                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>& VE_map,
                                             const shared_ptr<thread_pool> &pool) {
        while (!ED->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

            compute_delete_edge_set(G, ED, edge_truss_map, truss_order_map, rem, TS, VE_map, Ec, EU, pool);

            auto candidate_edge_map = make_shared<map<uint32_t,shared_ptr<vector<shared_ptr<abstract_edge>>>>>();
            for (const auto &[k,EUk]:*EU) {
                candidate_edge_map->insert({k, shared_ptr<vector<shared_ptr<abstract_edge>>>()});
            }

            /**
             * @brief parallel part
             */
            for (const auto &[k, EUk]:*EU) {
                pool->submit_task([=]() {
                    candidate_edge_map->at(k) = make_shared<vector<shared_ptr<abstract_edge>>>();
                    k_joint_delete(G, EUk, edge_truss_map, truss_order_map, rem, TS, VE_map, candidate_edge_map->at(k), k);
                });
            }
            pool->barrier();

            for(const auto&[k, k_vector]:*candidate_edge_map) {
                update_VE_map_for_removal(k_vector, edge_truss_map, VE_map);
            }

            for(const auto&[k, k_vector]:*candidate_edge_map) {
                pool->submit_task([=]{
                    for (const auto &e: *k_vector) {
                        truss_order_map->at(k)->remove(e);
                        --edge_truss_map->at(e);
                    }
                });
            }
            pool->barrier();

            for(const auto&[k, k_vector]:*candidate_edge_map) {
                update_VE_map_for_insertion(k_vector, edge_truss_map, VE_map);
            }

            for(const auto&[k, k_vector]:*candidate_edge_map) {
                pool->submit_task([=]{
                    for (const auto &e: *k_vector) {
                        truss_order_map->at(k - 1)->push_back(e);
                    }
                });
            }
            pool->barrier();

            for(const auto&[k, k_vector]:*candidate_edge_map)
            {
                pool->submit_task([=] {
                    for (const auto &e1: *k_vector) {
                        TS->at(e1) = 0;
                        rem->at(e1) = 0;

                        auto u = e1->get_source_vertex_id();
                        auto v = e1->get_destination_vertex_id();

                        if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                            swap(u, v);
                        }

                        auto u_map = VE_map->at(u);
                        for(auto iter = u_map->lower_bound(k - 1); iter!=u_map->end(); ++iter){
                            for(const auto &w:*iter->second){
                                auto e2 = G->get_edge(u, w);
                                if (w == v) {
                                    continue;
                                }
                                auto e3 = G->get_edge(v, w);
                                if (!e3 || edge_truss_map->at(e3) <  k - 1) {
                                    continue;
                                }

                                ++TS->at(e1);

                                if (test_order(edge_truss_map, truss_order_map, e1, e2) &&
                                    test_order(edge_truss_map, truss_order_map, e1, e3)) {
                                    ++rem->at(e1);
                                }
                            }
                        }

//                        for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
//                            if (w == v || edge_truss_map->at(e2) < k - 1) {
//                                continue;
//                            }
//                            auto e3 = G->get_edge(v, w);
//                            if (!e3 || edge_truss_map->at(e3) <  k - 1) {
//                                continue;
//                            }
//
//                            ++TS->at(e1);
//
//                            if (test_order(edge_truss_map, truss_order_map, e1, e2) &&
//                                test_order(edge_truss_map, truss_order_map, e1, e3)) {
//                                ++rem->at(e1);
//                            }
//                        }
                    }
                });
            }
            pool->barrier();
        }
    }


    void jes_order_truss_maintenance::k_joint_delete(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                                     uint32_t k) {
        auto TS = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e:*Ek){

            TS->insert({e, compute_truss_support(G, edge_truss_map, e, k)});
            if(TS->at(e) < k - 2){
                current_edge_set->insert(e);
            }
        }

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set){

                evicted_edge_set->insert(e1);
                candidate_edge_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u, v);
                }

                for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                    if(w == v || edge_truss_map->at(e2) < k || evicted_edge_set->count(e2)){
                        continue;
                    }
                    auto e3 = G->get_edge(v,w);
                    if(!e3 || edge_truss_map->at(e3) < k || evicted_edge_set->count(e3)){
                        continue;
                    }

                    if(edge_truss_map->at(e2) == k){
                        if(!TS->count(e2)){
                            TS->insert({e2, compute_truss_support(G, edge_truss_map, e2, k)});
                        }
                        --TS->at(e2);
                        if(TS->at(e2) < k - 2 && !current_edge_set->count(e2)){
                            next_edge_set->insert(e2);
                        }
                    }

                    if(edge_truss_map->at(e3) == k){
                        if(!TS->count(e3)){
                            TS->insert({e3, compute_truss_support(G, edge_truss_map, e3, k)});
                        }
                        --TS->at(e3);
                        if(TS->at(e3) < k - 2 && !current_edge_set->count(e3)){
                            next_edge_set->insert(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void jes_order_truss_maintenance::k_joint_delete(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                                     uint32_t k) {
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e:*Ek){
            if(TS->at(e) < k - 2){
                current_edge_set->insert(e);
            }
        }

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set){

                evicted_edge_set->insert(e1);
                candidate_edge_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u, v);
                }

                for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                    if(w == v || edge_truss_map->at(e2) < k || evicted_edge_set->count(e2)){
                        continue;
                    }
                    auto e3 = G->get_edge(v,w);
                    if(!e3 || edge_truss_map->at(e3) < k || evicted_edge_set->count(e3)){
                        continue;
                    }

                    if(edge_truss_map->at(e2) == k){
                        --TS->at(e2);
                        if(TS->at(e2) < k - 2 && !current_edge_set->count(e2)){
                            next_edge_set->insert(e2);
                        }
                    }

                    if(edge_truss_map->at(e3) == k){
                        --TS->at(e3);
                        if(TS->at(e3) < k - 2 && !current_edge_set->count(e3)){
                            next_edge_set->insert(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void jes_order_truss_maintenance::k_joint_delete(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                                     uint32_t k) {
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e:*Ek){
            if(TS->at(e) < k - 2){
                current_edge_set->insert(e);
            }
        }

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set){

                evicted_edge_set->insert(e1);
                candidate_edge_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u, v);
                }

                for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                    if(w == v || edge_truss_map->at(e2) < k || evicted_edge_set->count(e2)){
                        continue;
                    }
                    auto e3 = G->get_edge(v,w);
                    if(!e3 || edge_truss_map->at(e3) < k || evicted_edge_set->count(e3)){
                        continue;
                    }

                    if(edge_truss_map->at(e2) == k){
                        --TS->at(e2);
                        if(TS->at(e2) < k - 2 && !current_edge_set->count(e2)){
                            next_edge_set->insert(e2);
                        }
                    }

                    if(edge_truss_map->at(e3) == k){
                        --TS->at(e3);
                        if(TS->at(e3) < k - 2 && !current_edge_set->count(e3)){
                            next_edge_set->insert(e3);
                        }
                    }

                    auto min_e = e1;

                    if(test_order(edge_truss_map, truss_order_map, e2, min_e)){
                        min_e = e2;
                    }

                    if(test_order(edge_truss_map, truss_order_map, e3, min_e)){
                        min_e = e3;
                    }

                    if(min_e == e2){
                        --rem->at(e2);
                    }

                    if(min_e == e3){
                        --rem->at(e3);
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void jes_order_truss_maintenance::k_joint_delete(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<uint32_t,shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>& VE_map,
                                                     const shared_ptr<vector<shared_ptr<abstract_edge>>> &candidate_edge_vector,
                                                     uint32_t k) {
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e:*Ek){
            if(TS->at(e) < k - 2){
                current_edge_set->insert(e);
            }
        }

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set){

                evicted_edge_set->insert(e1);
                candidate_edge_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u, v);
                }

                auto u_map = VE_map->at(u);
                if(u_map->count(k)){
                    for(const auto &w:*u_map->at(k)){
                        auto e2 = G->get_edge(u, w);
                        if(w == v || evicted_edge_set->count(e2)){
                            continue;
                        }
                        auto e3 = G->get_edge(v,w);
                        if(!e3 || edge_truss_map->at(e3) < k || evicted_edge_set->count(e3)){
                            continue;
                        }

                        if(edge_truss_map->at(e2) == k){
                            --TS->at(e2);
                            if(TS->at(e2) < k - 2 && !current_edge_set->count(e2)){
                                next_edge_set->insert(e2);
                            }
                        }

                        if(edge_truss_map->at(e3) == k){
                            --TS->at(e3);
                            if(TS->at(e3) < k - 2 && !current_edge_set->count(e3)){
                                next_edge_set->insert(e3);
                            }
                        }

                        auto min_e = e1;

                        if(test_order(edge_truss_map, truss_order_map, e2, min_e)){
                            min_e = e2;
                        }

                        if(test_order(edge_truss_map, truss_order_map, e3, min_e)){
                            min_e = e3;
                        }

                        if(min_e == e2){
                            --rem->at(e2);
                        }

                        if(min_e == e3){
                            --rem->at(e3);
                        }
                    }
                }
                for(auto iter = u_map->lower_bound(k + 1); iter!=u_map->end(); ++iter){
                    for(const auto &w:*iter->second){
                        auto e2 = G->get_edge(u, w);
                        auto e3 = G->get_edge(v, w);
                        if(!e3 || edge_truss_map->at(e3) < k || evicted_edge_set->count(e3)){
                            continue;
                        }

                        if(edge_truss_map->at(e3) == k){
                            --TS->at(e3);
                            if(TS->at(e3) < k - 2 && !current_edge_set->count(e3)){
                                next_edge_set->insert(e3);
                            }
                        }

                        auto min_e = e1;

                        if(test_order(edge_truss_map, truss_order_map, e3, min_e)){
                            min_e = e3;
                        }

                        if(min_e == e3){
                            --rem->at(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }


    void jes_order_truss_maintenance::k_joint_insert(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                     uint32_t k) {
        auto current_edge_map = make_shared<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>>();
        auto current_order_list = truss_order_map->at(k);

        auto ext = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        auto e_edge_vector = make_shared<vector<shared_ptr<scnu::abstract_edge>>>();
        update_order_list(G, evicted_edge_set, candidate_edge_set, candidate_edge_support_map, edge_truss_map,
                          current_order_list, e_edge_vector, k);

        update_edge_support(G, candidate_edge_set, edge_truss_map, truss_order_map, ext, current_edge_map, k);

        for (auto iter = e_edge_vector->rbegin(); iter != e_edge_vector->rend(); ++iter) {
            current_order_list->left_insert(*iter);
        }

        while (!current_edge_map->empty()) {
            auto [key, p] = *current_edge_map->begin();
            current_edge_map->erase(key);

            auto p_next = p->get_next();
            auto e_current = p->get_value();

            auto rem_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto rem_value = get_rem_set(G, edge_truss_map, truss_order_map, e_current, rem_set, k);

            if (rem_value + ext->at(e_current) > k - 2) {
                /**
                 * @brief insert e_current into C_k
                 */
                candidate_edge_set->insert(e_current);

                candidate_edge_support_map->insert({e_current, ext->at(e_current) + rem_value});

                ext->at(e_current) = 0;

                for (const auto &e: *rem_set) {
                    if (current_order_list->count_value(e)) {
                        if(!ext->count(e)){
                            auto e_node = current_order_list->find(e);
                            current_edge_map->insert({e_node->get_key(), e_node});
                            ext->insert({e, 0});
                        }
                        ++ext->at(e);
                    }
                }

                /**
                 * @brief remove e_current from wing order
                 */
                current_order_list->remove(e_current);
            } else {
                ext->at(e_current) = 0;

                auto u = e_current->get_source_vertex_id();
                auto v = e_current->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || !(candidate_edge_set->count(e2) || edge_truss_map->at(e2) > k || current_order_list->count_value(e2) &&
                                                                                                           test_order(current_order_list, e_current, e2))) {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if (!e3 || !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k || current_order_list->count_value(e3) &&
                                                                                                test_order(current_order_list, e_current, e3))) {
                        continue;
                    }

                    if (!candidate_edge_set->count(e2) && !candidate_edge_set->count(e3)) {
                        continue;
                    }

                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2) {
                            evicted_edge_set->insert(e2);
                        }
                    } else if (ext->count(e2)) {
                        --ext->at(e2);
                        if (ext->at(e2) == 0) {
                            ext->erase(e2);
                            auto e2_node = current_order_list->find(e2);
                            current_edge_map->erase(e2_node->get_key());
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2) {
                            evicted_edge_set->insert(e3);
                        }
                    } else if (ext->count(e3)) {
                        --ext->at(e3);
                        if (ext->at(e3) == 0) {
                            ext->erase(e3);
                            auto e3_node = current_order_list->find(e3);
                            current_edge_map->erase(e3_node->get_key());
                        }
                    }
                }
                remove_unsatisfied_edges(G, evicted_edge_set, candidate_edge_set, edge_truss_map, truss_order_map,
                                         candidate_edge_support_map, ext, current_edge_map,
                                         e_current, k);
            }
            p = p_next;
        }

        current_order_list->reset_order();
    }

    void jes_order_truss_maintenance::k_joint_insert(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     uint32_t k) {
        auto current_edge_map = make_shared<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>>();
        auto current_order_list = truss_order_map->at(k);

        auto ext = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        auto e_edge_vector = make_shared<vector<shared_ptr<scnu::abstract_edge>>>();
        update_order_list(G, evicted_edge_set, candidate_edge_set, candidate_edge_support_map, edge_truss_map,
                          current_order_list, TS, e_edge_vector, k);

        update_edge_support(G, candidate_edge_set, edge_truss_map, truss_order_map, ext, current_edge_map, k);

        for (auto iter = e_edge_vector->rbegin(); iter != e_edge_vector->rend(); ++iter) {
            current_order_list->left_insert(*iter);
        }

        while (!current_edge_map->empty()) {
            auto [key, p] = *current_edge_map->begin();
            current_edge_map->erase(key);

            auto p_next = p->get_next();
            auto e_current = p->get_value();

            auto rem_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto rem_value = get_rem_set(G, edge_truss_map, truss_order_map, e_current, rem_set, k);

            if (rem_value + ext->at(e_current) > k - 2) {
                /**
                 * @brief insert e_current into C_k
                 */
                candidate_edge_set->insert(e_current);

                candidate_edge_support_map->insert({e_current, ext->at(e_current) + rem_value});

                ext->at(e_current) = 0;

                for (const auto &e: *rem_set) {
                    if (current_order_list->count_value(e)) {
                        if(!ext->count(e)){
                            auto e_node = current_order_list->find(e);
                            current_edge_map->insert({e_node->get_key(), e_node});
                            ext->insert({e, 0});
                        }
                        ++ext->at(e);
                    }
                }

                /**
                 * @brief remove e_current from wing order
                 */
                current_order_list->remove(e_current);
            } else {
                ext->at(e_current) = 0;

                auto u = e_current->get_source_vertex_id();
                auto v = e_current->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || !(candidate_edge_set->count(e2) || edge_truss_map->at(e2) > k || current_order_list->count_value(e2) &&
                                                                                                   test_order(current_order_list, e_current, e2))) {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if (!e3 || !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k || current_order_list->count_value(e3) &&
                                                                                                test_order(current_order_list, e_current, e3))) {
                        continue;
                    }

                    if (!candidate_edge_set->count(e2) && !candidate_edge_set->count(e3)) {
                        continue;
                    }

                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2) {
                            evicted_edge_set->insert(e2);
                        }
                    } else if (ext->count(e2)) {
                        --ext->at(e2);
                        if (ext->at(e2) == 0) {
                            ext->erase(e2);
                            auto e2_node = current_order_list->find(e2);
                            current_edge_map->erase(e2_node->get_key());
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2) {
                            evicted_edge_set->insert(e3);
                        }
                    } else if (ext->count(e3)) {
                        --ext->at(e3);
                        if (ext->at(e3) == 0) {
                            ext->erase(e3);
                            auto e3_node = current_order_list->find(e3);
                            current_edge_map->erase(e3_node->get_key());
                        }
                    }
                }
                remove_unsatisfied_edges(G, evicted_edge_set, candidate_edge_set, edge_truss_map, truss_order_map,
                                         candidate_edge_support_map, ext, current_edge_map,
                                         e_current, k);
            }
            p = p_next;
        }

        current_order_list->reset_order();
    }

    void jes_order_truss_maintenance::k_joint_insert(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     uint32_t k) {
        auto current_edge_map = make_shared<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>>();
        auto current_order_list = truss_order_map->at(k);

        auto ext = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        auto e_edge_vector = make_shared<vector<shared_ptr<scnu::abstract_edge>>>();

        update_order_list(G, evicted_edge_set, candidate_edge_set, candidate_edge_support_map, edge_truss_map,
                          current_order_list, rem, TS, e_edge_vector, k);

        update_edge_support(G, candidate_edge_set, edge_truss_map, truss_order_map, ext, current_edge_map, k);

        for (auto iter = e_edge_vector->rbegin(); iter != e_edge_vector->rend(); ++iter) {
            current_order_list->left_insert(*iter);
        }
//
//        for(auto p = current_order_list->get_head(); p;p=p->get_next()){
//            auto e = p->get_value();
//            printf("(%u,%u) ", e->get_source_vertex_id(), e->get_destination_vertex_id());
//        }
//        printf("\n");

        while (!current_edge_map->empty()) {
            auto [key, p] = *current_edge_map->begin();
            current_edge_map->erase(key);

            auto p_next = p->get_next();
            auto e_current = p->get_value();


            if (rem->at(e_current) + ext->at(e_current) > k - 2) {
                /**
                 * @brief insert e_current into C_k
                 */
                candidate_edge_set->insert(e_current);
                auto rem_set  = get_rem_set(G, edge_truss_map, current_order_list, e_current, k);

                candidate_edge_support_map->insert({e_current, ext->at(e_current) + rem->at(e_current)});

                ext->at(e_current) = 0;

                for (const auto &e: *rem_set) {
                    if (current_order_list->count_value(e)) {
                        if(!ext->count(e)){
                            auto e_node = current_order_list->find(e);
                            current_edge_map->insert({e_node->get_key(), e_node});
                            ext->insert({e, 0});
                        }
                        ++ext->at(e);
                    }
                }

                /**
                 * @brief remove e_current from wing order
                 */
                current_order_list->remove(e_current);
            } else {
                rem->at(e_current) = rem->at(e_current) + ext->at(e_current);
                ext->at(e_current) = 0;

                auto u = e_current->get_source_vertex_id();
                auto v = e_current->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k ||
                        !(candidate_edge_set->count(e2) || edge_truss_map->at(e2) > k ||
                          current_order_list->count_value(e2) &&
                          test_order(current_order_list, e_current, e2))) {
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k ||
                        !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                          current_order_list->count_value(e3) &&
                          test_order(current_order_list, e_current, e3))) {
                        continue;
                    }

                    if (!candidate_edge_set->count(e2) && !candidate_edge_set->count(e3)) {
                        continue;
                    }

                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2) {
                            evicted_edge_set->insert(e2);
                        }
                    } else if (ext->count(e2)) {
                        --ext->at(e2);
                        if (ext->at(e2) == 0) {
                            ext->erase(e2);
                            auto e2_node = current_order_list->find(e2);
                            current_edge_map->erase(e2_node->get_key());
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2) {
                            evicted_edge_set->insert(e3);
                        }
                    } else if (ext->count(e3)) {
                        --ext->at(e3);
                        if (ext->at(e3) == 0) {
                            ext->erase(e3);
                            auto e3_node = current_order_list->find(e3);
                            current_edge_map->erase(e3_node->get_key());
                        }
                    }
                }
                remove_unsatisfied_edges(G, evicted_edge_set, candidate_edge_set, edge_truss_map, truss_order_map, rem,
                                         candidate_edge_support_map, ext, current_edge_map,
                                         e_current, k);
            }
            p = p_next;
        }

        current_order_list->reset_order();
    }

    void jes_order_truss_maintenance::k_joint_insert(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                     const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                                     uint32_t k) {
        auto current_edge_map = make_shared<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>>();
        auto current_order_list = truss_order_map->at(k);

        auto ext = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        auto e_edge_vector = make_shared<vector<shared_ptr<scnu::abstract_edge>>>();

        update_order_list(G, evicted_edge_set, candidate_edge_set, candidate_edge_support_map, edge_truss_map,
                          current_order_list, rem, TS, VE_map, e_edge_vector, k);

        update_edge_support(G, candidate_edge_set, edge_truss_map, truss_order_map, VE_map, ext, current_edge_map, k);

        for (auto iter = e_edge_vector->rbegin(); iter != e_edge_vector->rend(); ++iter) {
            current_order_list->left_insert(*iter);
        }
//
//        for(auto p = current_order_list->get_head(); p;p=p->get_next()){
//            auto e = p->get_value();
//            printf("(%u,%u) ", e->get_source_vertex_id(), e->get_destination_vertex_id());
//        }
//        printf("\n");

        while (!current_edge_map->empty()) {
            auto [key, p] = *current_edge_map->begin();
            current_edge_map->erase(key);

            auto p_next = p->get_next();
            auto e_current = p->get_value();


            if (rem->at(e_current) + ext->at(e_current) > k - 2) {
                /**
                 * @brief insert e_current into C_k
                 */
                candidate_edge_set->insert(e_current);

                auto rem_set = get_rem_set(G, edge_truss_map, current_order_list, VE_map, e_current, k);

                candidate_edge_support_map->insert({e_current, ext->at(e_current) + rem->at(e_current)});

                ext->at(e_current) = 0;

                for (const auto &e: *rem_set) {
                    if (current_order_list->count_value(e)) {
                        if (!ext->count(e)) {
                            auto e_node = current_order_list->find(e);
                            current_edge_map->insert({e_node->get_key(), e_node});
                            ext->insert({e, 0});
                        }
                        ++ext->at(e);
                    }
                }

                /**
                 * @brief remove e_current from wing order
                 */
                current_order_list->remove(e_current);
            } else {
                rem->at(e_current) = rem->at(e_current) + ext->at(e_current);
                ext->at(e_current) = 0;

                auto u = e_current->get_source_vertex_id();
                auto v = e_current->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                auto u_map = VE_map->at(u);
                if(u_map->count(k)){
                    for(const auto&w:*u_map->at(k)){
                        auto e2 = G->get_edge(u, w);
                        if (w == v || !(candidate_edge_set->count(e2) ||
                                        current_order_list->count_value(e2) &&
                                        test_order(current_order_list, e_current, e2))) {
                            continue;
                        }
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k ||
                            !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                              current_order_list->count_value(e3) &&
                              test_order(current_order_list, e_current, e3))) {
                            continue;
                        }

                        if (!candidate_edge_set->count(e2) && !candidate_edge_set->count(e3)) {
                            continue;
                        }

                        if (candidate_edge_set->count(e2)) {
                            --candidate_edge_support_map->at(e2);
                            if (candidate_edge_support_map->at(e2) <= k - 2) {
                                evicted_edge_set->insert(e2);
                            }
                        } else if (ext->count(e2)) {
                            --ext->at(e2);
                            if (ext->at(e2) == 0) {
                                ext->erase(e2);
                                auto e2_node = current_order_list->find(e2);
                                current_edge_map->erase(e2_node->get_key());
                            }
                        }

                        if (candidate_edge_set->count(e3)) {
                            --candidate_edge_support_map->at(e3);
                            if (candidate_edge_support_map->at(e3) <= k - 2) {
                                evicted_edge_set->insert(e3);
                            }
                        } else if (ext->count(e3)) {
                            --ext->at(e3);
                            if (ext->at(e3) == 0) {
                                ext->erase(e3);
                                auto e3_node = current_order_list->find(e3);
                                current_edge_map->erase(e3_node->get_key());
                            }
                        }
                    }
                }
                for(auto iter = u_map->lower_bound(k + 1); iter != u_map->end(); ++iter){
                    for(const auto&w:*iter->second){
                        auto e2 = G->get_edge(u, w);
                        auto e3 = G->get_edge(v, w);

                        if (!e3 || edge_truss_map->at(e3) < k ||
                            !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                              current_order_list->count_value(e3) &&
                              test_order(current_order_list, e_current, e3))) {
                            continue;
                        }

                        if (!candidate_edge_set->count(e2) && !candidate_edge_set->count(e3)) {
                            continue;
                        }

                        if (candidate_edge_set->count(e3)) {
                            --candidate_edge_support_map->at(e3);
                            if (candidate_edge_support_map->at(e3) <= k - 2) {
                                evicted_edge_set->insert(e3);
                            }
                        } else if (ext->count(e3)) {
                            --ext->at(e3);
                            if (ext->at(e3) == 0) {
                                ext->erase(e3);
                                auto e3_node = current_order_list->find(e3);
                                current_edge_map->erase(e3_node->get_key());
                            }
                        }
                    }
                }

//                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
//                    if (w == v || edge_truss_map->at(e2) < k ||
//                        !(candidate_edge_set->count(e2) || edge_truss_map->at(e2) > k ||
//                          current_order_list->count_value(e2) &&
//                          test_order(current_order_list, e_current, e2))) {
//                        continue;
//                    }
//                    auto e3 = G->get_edge(v, w);
//                    if (!e3 || edge_truss_map->at(e3) < k ||
//                        !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
//                          current_order_list->count_value(e3) &&
//                          test_order(current_order_list, e_current, e3))) {
//                        continue;
//                    }
//
//                    if (!candidate_edge_set->count(e2) && !candidate_edge_set->count(e3)) {
//                        continue;
//                    }
//
//                    if (candidate_edge_set->count(e2)) {
//                        --candidate_edge_support_map->at(e2);
//                        if (candidate_edge_support_map->at(e2) <= k - 2) {
//                            evicted_edge_set->insert(e2);
//                        }
//                    } else if (ext->count(e2)) {
//                        --ext->at(e2);
//                        if (ext->at(e2) == 0) {
//                            ext->erase(e2);
//                            auto e2_node = current_order_list->find(e2);
//                            current_edge_map->erase(e2_node->get_key());
//                        }
//                    }
//
//                    if (candidate_edge_set->count(e3)) {
//                        --candidate_edge_support_map->at(e3);
//                        if (candidate_edge_support_map->at(e3) <= k - 2) {
//                            evicted_edge_set->insert(e3);
//                        }
//                    } else if (ext->count(e3)) {
//                        --ext->at(e3);
//                        if (ext->at(e3) == 0) {
//                            ext->erase(e3);
//                            auto e3_node = current_order_list->find(e3);
//                            current_edge_map->erase(e3_node->get_key());
//                        }
//                    }
//                }
                remove_unsatisfied_edges(G, evicted_edge_set, candidate_edge_set, edge_truss_map, truss_order_map, rem, VE_map,
                                         candidate_edge_support_map, ext, current_edge_map,
                                         e_current, k);
            }
            p = p_next;
        }

        current_order_list->reset_order();
    }

    void jes_order_truss_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                                               const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                                               const shared_ptr<abstract_edge> &e_pivot,
                                                               uint32_t k) {
        auto current_order_list = truss_order_map->at(k);
        auto e_vector = make_shared<vector<shared_ptr<abstract_edge>>>();

        auto e_pivot_node = current_order_list->find(e_pivot)->get_next();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *current_edge_set) {
                candidate_edge_set->erase(e1);

                e_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k ||
                        !(candidate_edge_set->count(e2) || edge_truss_map->at(e2) > k ||
                          current_order_list->count_value(e2) &&
                          test_order(current_order_list, e_pivot, e2))) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k ||
                        !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                          current_order_list->count_value(e3) &&
                          test_order(current_order_list, e_pivot, e3))) {
                        continue;
                    }

                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                            next_edge_set->insert(e2);
                        }
                    } else if (ext->count(e2)) {
                        --ext->at(e2);
                        if (ext->at(e2) == 0) {
                            ext->erase(e2);
                            auto e2_node = current_order_list->find(e2);
                            current_edge_map->erase(e2_node->get_key());
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                            next_edge_set->insert(e3);
                        }
                    } else if (ext->count(e3)) {
                        --ext->at(e3);
                        if (ext->at(e3) == 0) {
                            ext->erase(e3);
                            auto e3_node = current_order_list->find(e3);
                            current_edge_map->erase(e3_node->get_key());
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*e_vector){
            auto e_node = make_shared<extend_node<int, shared_ptr<abstract_edge>>>(0, e);
            if (!e_pivot_node) {
                current_order_list->right_insert(e_node);
            } else {
                current_order_list->insert_before(e_node, e_pivot_node);
            }
        }
    }

    void jes_order_truss_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                                               const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                                               const shared_ptr<abstract_edge> &e_pivot,
                                                               uint32_t k) {
        auto current_order_list = truss_order_map->at(k);
        auto e_vector = make_shared<vector<shared_ptr<abstract_edge>>>();

        auto e_pivot_node = current_order_list->find(e_pivot)->get_next();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *current_edge_set) {
                candidate_edge_set->erase(e1);

                e_vector->push_back(e1);
                rem->at(e1) =  candidate_edge_support_map->at(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k ||
                        !(candidate_edge_set->count(e2) || edge_truss_map->at(e2) > k ||
                          current_order_list->count_value(e2) &&
                          test_order(current_order_list, e_pivot, e2))) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k ||
                        !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                          current_order_list->count_value(e3) &&
                          test_order(current_order_list, e_pivot, e3))) {
                        continue;
                    }

                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                            next_edge_set->insert(e2);
                        }
                    } else if (ext->count(e2)) {
                        --ext->at(e2);
                        if (ext->at(e2) == 0) {
                            ext->erase(e2);
                            auto e2_node = current_order_list->find(e2);
                            current_edge_map->erase(e2_node->get_key());
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                            next_edge_set->insert(e3);
                        }
                    } else if (ext->count(e3)) {
                        --ext->at(e3);
                        if (ext->at(e3) == 0) {
                            ext->erase(e3);
                            auto e3_node = current_order_list->find(e3);
                            current_edge_map->erase(e3_node->get_key());
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*e_vector){
            auto e_node = make_shared<extend_node<int, shared_ptr<abstract_edge>>>(0, e);
            if (!e_pivot_node) {
                current_order_list->right_insert(e_node);
            } else {
                current_order_list->insert_before(e_node, e_pivot_node);
            }
        }
    }

    void jes_order_truss_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>& truss_order_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                               const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> & VE_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                                               const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                                               const shared_ptr<abstract_edge> &e_pivot,
                                                               uint32_t k) {
        auto current_order_list = truss_order_map->at(k);
        auto e_vector = make_shared<vector<shared_ptr<abstract_edge>>>();

        auto e_pivot_node = current_order_list->find(e_pivot)->get_next();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *current_edge_set) {
                candidate_edge_set->erase(e1);

                e_vector->push_back(e1);
                rem->at(e1) =  candidate_edge_support_map->at(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                auto u_map = VE_map->at(u);
                if(u_map->count(k)){
                    for(const auto &w:*u_map->at(k)){
                        auto e2 = G->get_edge(u, w);
                        if (w == v ||
                            !(candidate_edge_set->count(e2) ||
                              current_order_list->count_value(e2) &&
                              test_order(current_order_list, e_pivot, e2))) {
                            continue;
                        }

                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k ||
                            !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                              current_order_list->count_value(e3) &&
                              test_order(current_order_list, e_pivot, e3))) {
                            continue;
                        }

                        if (candidate_edge_set->count(e2)) {
                            --candidate_edge_support_map->at(e2);
                            if (candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                                next_edge_set->insert(e2);
                            }
                        } else if (ext->count(e2)) {
                            --ext->at(e2);
                            if (ext->at(e2) == 0) {
                                ext->erase(e2);
                                auto e2_node = current_order_list->find(e2);
                                current_edge_map->erase(e2_node->get_key());
                            }
                        }

                        if (candidate_edge_set->count(e3)) {
                            --candidate_edge_support_map->at(e3);
                            if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                                next_edge_set->insert(e3);
                            }
                        } else if (ext->count(e3)) {
                            --ext->at(e3);
                            if (ext->at(e3) == 0) {
                                ext->erase(e3);
                                auto e3_node = current_order_list->find(e3);
                                current_edge_map->erase(e3_node->get_key());
                            }
                        }
                    }
                }
                for(auto iter = u_map->lower_bound(k + 1); iter!=u_map->end(); ++iter){
                    for(const auto &w:*iter->second){
                        auto e2 = G->get_edge(u, w);
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k ||
                            !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                              current_order_list->count_value(e3) &&
                              test_order(current_order_list, e_pivot, e3))) {
                            continue;
                        }

                        if (candidate_edge_set->count(e3)) {
                            --candidate_edge_support_map->at(e3);
                            if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                                next_edge_set->insert(e3);
                            }
                        } else if (ext->count(e3)) {
                            --ext->at(e3);
                            if (ext->at(e3) == 0) {
                                ext->erase(e3);
                                auto e3_node = current_order_list->find(e3);
                                current_edge_map->erase(e3_node->get_key());
                            }
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }

        for(const auto &e:*e_vector){
            auto e_node = make_shared<extend_node<int, shared_ptr<abstract_edge>>>(0, e);
            if (!e_pivot_node) {
                current_order_list->right_insert(e_node);
            } else {
                current_order_list->insert_before(e_node, e_pivot_node);
            }
        }
    }

    bool jes_order_truss_maintenance::test_order(
            const shared_ptr<unordered_map<shared_ptr<scnu::abstract_edge>, uint32_t>> &edge_truss_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>>&truss_order_map,
            const shared_ptr<abstract_edge> &e1, const shared_ptr<scnu::abstract_edge> &e2) {
        if(edge_truss_map->at(e1) < edge_truss_map->at(e2)){
            return true;
        }else if(edge_truss_map->at(e1) > edge_truss_map->at(e2)){
            return false;
        }else{
            auto order_list = truss_order_map->at(edge_truss_map->at(e1));
            return order_list->find(e1)->get_key() < order_list->find(e2)->get_key();
        }
    }

    bool jes_order_truss_maintenance::test_order(
            const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>&order_list,
            const shared_ptr<abstract_edge> &e1, const shared_ptr<scnu::abstract_edge> &e2) {
        return order_list->find(e1)->get_key() < order_list->find(e2)->get_key();
    }

    void jes_order_truss_maintenance::update_edge_support(const shared_ptr<abstract_graph> &G,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                                          const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                                          uint32_t k) {
        auto current_order_list = truss_order_map->at(k);

        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        for (const auto &e1: *candidate_edge_set) {
            visited_edge_set->insert(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || visited_edge_set->count(e2) || edge_truss_map->at(e2) < k ||
                    !(candidate_edge_set->count(e2) || edge_truss_map->at(e2) > k ||
                      current_order_list->count_value(e2))) {
                    continue;
                }


                auto e3 = G->get_edge(v, w);
                if (!e3 || visited_edge_set->count(e3) || edge_truss_map->at(e3) < k ||
                    !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                      current_order_list->count_value(e3))) {
                    continue;
                }


                if (current_order_list->count_value(e2)) {
                    if(!ext->count(e2)){
                        ext->insert({e2, 0});
                        auto e2_node = current_order_list->find(e2);
                        current_edge_map->insert({e2_node->get_key(), e2_node});
                    }
                    ++ext->at(e2);
                }

                if (current_order_list->count_value(e3)) {
                    if(!ext->count(e3)){
                        ext->insert({e3, 0});
                        auto e3_node = current_order_list->find(e3);
                        current_edge_map->insert({e3_node->get_key(), e3_node});
                    }
                    ++ext->at(e3);
                }
            }
        }
    }

    void jes_order_truss_maintenance::update_edge_support(const shared_ptr<abstract_graph> &G,
                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<extend_list<int, shared_ptr<abstract_edge>>>>> &truss_order_map,
                                                          const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>>& VE_map,
                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &ext,
                                                          const shared_ptr<map<int, shared_ptr<extend_node<int, shared_ptr<abstract_edge>>>>> &current_edge_map,
                                                          uint32_t k) {
        auto current_order_list = truss_order_map->at(k);

        auto visited_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

        for (const auto &e1: *candidate_edge_set) {
            visited_edge_set->insert(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            auto u_map = VE_map->at(u);
            if(u_map->count(k)){
                for(const auto &w:*u_map->at(k)){
                    auto e2 = G->get_edge(u, w);
                    if(w == v || visited_edge_set->count(e2) ||  !(candidate_edge_set->count(e2) ||
                                                                   current_order_list->count_value(e2))){
                        continue;
                    }
                    auto e3 = G->get_edge(v, w);

                    if (!e3 || visited_edge_set->count(e3) || edge_truss_map->at(e3) < k ||
                        !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                          current_order_list->count_value(e3))) {
                        continue;
                    }

                    if (current_order_list->count_value(e2)) {
                        if(!ext->count(e2)){
                            ext->insert({e2, 0});
                            auto e2_node = current_order_list->find(e2);
                            current_edge_map->insert({e2_node->get_key(), e2_node});
                        }
                        ++ext->at(e2);
                    }

                    if (current_order_list->count_value(e3)) {
                        if(!ext->count(e3)){
                            ext->insert({e3, 0});
                            auto e3_node = current_order_list->find(e3);
                            current_edge_map->insert({e3_node->get_key(), e3_node});
                        }
                        ++ext->at(e3);
                    }
                }
            }

            for(auto iter = u_map->lower_bound(k + 1); iter!=u_map->end(); ++iter){
                for(const auto &w:*iter->second){
                    auto e2 = G->get_edge(u, w);
                    auto e3 = G->get_edge(v, w);

                    if (!e3 || visited_edge_set->count(e3) || edge_truss_map->at(e3) < k ||
                        !(candidate_edge_set->count(e3) || edge_truss_map->at(e3) > k ||
                          current_order_list->count_value(e3))) {
                        continue;
                    }


                    if (current_order_list->count_value(e3)) {
                        if(!ext->count(e3)){
                            ext->insert({e3, 0});
                            auto e3_node = current_order_list->find(e3);
                            current_edge_map->insert({e3_node->get_key(), e3_node});
                        }
                        ++ext->at(e3);
                    }
                }
            }
        }
    }

    void jes_order_truss_maintenance::update_edge_truss_support(const shared_ptr<abstract_graph> &G,
                                                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                                uint32_t k) {
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();


        for (const auto &e1: *edge_set) {
            visited_set->insert(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for(const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v || visited_set->count(e2) || edge_truss_map->at(e2) < k){
                    continue;
                }

                auto e3 = G->get_edge(v,w);
                if(!e3 || visited_set->count(e3) || edge_truss_map->at(e3) < k){
                    continue;
                }

                if (!edge_set->count(e2) && edge_truss_map->at(e2) == k) {
                    ++TS->at(e2);
                }

                if (!edge_set->count(e3) && edge_truss_map->at(e3) == k) {
                    ++TS->at(e3);
                }
            }
        }
    }

    void jes_order_truss_maintenance::update_edge_truss_support(const shared_ptr<abstract_graph> &G,
                                                                const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                                                uint32_t k) {
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();


        for (const auto &e1: *edge_set) {
            visited_set->insert(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            auto u_map = VE_map->at(u);
            for(auto iter = u_map->lower_bound(k); iter != u_map->end(); ++iter){
                for(const auto &w:*iter->second){
                    auto e2 = G->get_edge(u, w);
                    if(w == v || visited_set->count(e2)){
                        continue;
                    }

                    auto e3 = G->get_edge(v,w);
                    if(!e3 || visited_set->count(e3) || edge_truss_map->at(e3) < k){
                        continue;
                    }

                    if (!edge_set->count(e2) && edge_truss_map->at(e2) == k) {
                        ++TS->at(e2);
                    }

                    if (!edge_set->count(e3) && edge_truss_map->at(e3) == k) {
                        ++TS->at(e3);
                    }
                }
            }
        }
    }

    void jes_order_truss_maintenance::update_order_list(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                        const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                                        const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                                        uint32_t k) {
        auto e0_node = current_order_list->get_head();
        for (const auto &e1: *candidate_edge_set) {
            candidate_edge_support_map->insert({e1, 0});
            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                swap(u, v);
            }

            for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                if (w == v || edge_truss_map->at(e2) < k) {
                    continue;
                }
                auto e3 = G->get_edge(v, w);
                if (!e3 || edge_truss_map->at(e3) < k) {
                    continue;
                }

                ++candidate_edge_support_map->at(e1);
            }

            if (candidate_edge_support_map->at(e1) <= k - 2) {
                current_edge_set->insert(e1);
            }
        }

        auto evicted_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *current_edge_set) {
                candidate_edge_set->erase(e1);
                evicted_set->insert(e1);

                e_edge_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k || evicted_set->count(e2)) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k || evicted_set->count(e3)) {
                        continue;
                    }


                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                            next_edge_set->insert(e2);
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                            next_edge_set->insert(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void jes_order_truss_maintenance::update_order_list(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                        const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                        const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                                        uint32_t k) {
        auto e0_node = current_order_list->get_head();
        for (const auto &e1: *candidate_edge_set) {
            candidate_edge_support_map->insert({e1, TS->at(e1)});

            if (candidate_edge_support_map->at(e1) <= k - 2) {
                current_edge_set->insert(e1);
            }
        }

        auto evicted_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *current_edge_set) {
                candidate_edge_set->erase(e1);
                evicted_set->insert(e1);

                e_edge_vector->push_back(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k || evicted_set->count(e2)) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k || evicted_set->count(e3)) {
                        continue;
                    }


                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                            next_edge_set->insert(e2);
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                            next_edge_set->insert(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void jes_order_truss_maintenance::update_order_list(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                        const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                        const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                                        uint32_t k) {
        auto e0_node = current_order_list->get_head();
        for (const auto &e1: *candidate_edge_set) {
            candidate_edge_support_map->insert({e1, TS->at(e1)});

            if (candidate_edge_support_map->at(e1) <= k - 2) {
                current_edge_set->insert(e1);
            }
        }

        auto evicted_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *current_edge_set) {
                candidate_edge_set->erase(e1);
                evicted_set->insert(e1);

                e_edge_vector->push_back(e1);
                rem->at(e1) = candidate_edge_support_map->at(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                for (const auto &[w, e2]: *G->get_vertex(u)->get_edge_map()) {
                    if (w == v || edge_truss_map->at(e2) < k || evicted_set->count(e2)) {
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if (!e3 || edge_truss_map->at(e3) < k || evicted_set->count(e3)) {
                        continue;
                    }


                    if (candidate_edge_set->count(e2)) {
                        --candidate_edge_support_map->at(e2);
                        if (candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                            next_edge_set->insert(e2);
                        }
                    }

                    if (candidate_edge_set->count(e3)) {
                        --candidate_edge_support_map->at(e3);
                        if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                            next_edge_set->insert(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void jes_order_truss_maintenance::update_order_list(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                        const shared_ptr<extend_list<int, shared_ptr<abstract_edge>>> &current_order_list,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &rem,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &TS,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>>> &VE_map,
                                                        const shared_ptr<vector<shared_ptr<scnu::abstract_edge>>> &e_edge_vector,
                                                        uint32_t k) {
        auto e0_node = current_order_list->get_head();
        for (const auto &e1: *candidate_edge_set) {
            candidate_edge_support_map->insert({e1, TS->at(e1)});

            if (candidate_edge_support_map->at(e1) <= k - 2) {
                current_edge_set->insert(e1);
            }
        }

        auto evicted_set = make_shared<unordered_set<shared_ptr<scnu::abstract_edge>>>();
        while (!current_edge_set->empty()) {
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for (const auto &e1: *current_edge_set) {
                candidate_edge_set->erase(e1);
                evicted_set->insert(e1);

                e_edge_vector->push_back(e1);
                rem->at(e1) = candidate_edge_support_map->at(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if (G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()) {
                    swap(u, v);
                }

                auto u_map = VE_map->at(u);
                if(u_map->count(k)){
                    for (const auto &w: *u_map->at(k)) {
                        auto e2 = G->get_edge(u, w);
                        if (w == v || evicted_set->count(e2)) {
                            continue;
                        }
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k || evicted_set->count(e3)) {
                            continue;
                        }

                        if (candidate_edge_set->count(e2)) {
                            --candidate_edge_support_map->at(e2);
                            if (candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)) {
                                next_edge_set->insert(e2);
                            }
                        }

                        if (candidate_edge_set->count(e3)) {
                            --candidate_edge_support_map->at(e3);
                            if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                                next_edge_set->insert(e3);
                            }
                        }
                    }
                }
                for(auto iter = u_map->lower_bound(k + 1); iter != u_map->end(); ++iter){
                    for(const auto&w:*iter->second){
                        auto e2 = G->get_edge(u, w);
                        auto e3 = G->get_edge(v, w);
                        if (!e3 || edge_truss_map->at(e3) < k || evicted_set->count(e3)) {
                            continue;
                        }

                        if (candidate_edge_set->count(e3)) {
                            --candidate_edge_support_map->at(e3);
                            if (candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)) {
                                next_edge_set->insert(e3);
                            }
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }
}