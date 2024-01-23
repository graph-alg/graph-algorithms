
#include "truss/jes_truss_maintenance.h"

namespace scnu {
    uint32_t jes_truss_maintenance::compute_truss_number(const shared_ptr<abstract_graph> &G,
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

            uint32_t min_k = std::min(edge_truss_map->at(e2),edge_truss_map->at(e3));
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

    void jes_truss_maintenance::compute_delete_edge_set(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU) {

        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto superior_edge_triangle_map = make_shared<unordered_map<shared_ptr<abstract_edge>, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        compute_delete_superior_edge_set(G, ED, edge_truss_map, superior_edge_set, superior_edge_triangle_map);
        compute_maximal_three_hop_independent_set(G, superior_edge_set, superior_edge_triangle_map, Ec);
        /**
         * @brief parallel part
         */
        auto mark = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
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

                if (edge_truss_map->at(e2) <= edge_truss_map->at(e1)) {
                    if (!EU->count(edge_truss_map->at(e2))) {
                        EU->insert({edge_truss_map->at(e2), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                    }
                    EU->at(edge_truss_map->at(e2))->insert(e2);
                }

                if (edge_truss_map->at(e3) <= edge_truss_map->at(e1)) {
                    if (!EU->count(edge_truss_map->at(e3))) {
                        EU->insert({edge_truss_map->at(e3), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                    }
                    EU->at(edge_truss_map->at(e3))->insert(e3);
                }
            }
        }

        for(const auto &e:*Ec){
            G->remove_edge(e);
            ED->erase(e);
            edge_truss_map->erase(e);
        }
    }

    void jes_truss_maintenance::compute_insert_edge_set(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &EU) {
        /**
          * @brief get G+EI
          * @remarks  next loop will operate on current graph
        */
        G->insert_edge_collection(EI);
        auto superior_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto superior_edge_triangle_map = make_shared<unordered_map<shared_ptr<abstract_edge>,
                shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();

        compute_insert_superior_edge_set(G, EI, edge_truss_map, superior_edge_set, superior_edge_triangle_map);
        compute_maximal_three_hop_independent_set(G, superior_edge_set,superior_edge_triangle_map, Ec);
        G->remove_edge_collection(EI);

        for(const auto &e:*Ec){
            G->insert_edge(e);
            EI->erase(e);
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

                if (edge_truss_map->at(e2) <= pre_truss_number) {
                    if (!EU->count(edge_truss_map->at(e2))) {
                        EU->insert({edge_truss_map->at(e2), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                    }
                    EU->at(edge_truss_map->at(e2))->insert(e2);

                    if(edge_truss_map->at(e2) == pre_truss_number){
                        EU->at(edge_truss_map->at(e2))->insert(e1);
                    }
                }

                if (edge_truss_map->at(e3) <= pre_truss_number) {
                    if (!EU->count(edge_truss_map->at(e3))) {
                        EU->insert({edge_truss_map->at(e3), make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
                    }
                    EU->at(edge_truss_map->at(e3))->insert(e3);

                    if(edge_truss_map->at(e3) == pre_truss_number){
                        EU->at(edge_truss_map->at(e3))->insert(e1);
                    }
                }
            }
        }
    }

    uint32_t
    jes_truss_maintenance::compute_pre_truss_number(const shared_ptr<abstract_graph> &G,
                                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                    const shared_ptr<abstract_edge>& e) {
        return compute_truss_number(G, edge_truss_map, e);
    }


    void jes_truss_maintenance::compute_maximal_three_hop_independent_set(const shared_ptr<abstract_graph> &G,
                                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &superior_edge_set,
                                                                          const shared_ptr<unordered_map<shared_ptr<abstract_edge>,
                                                                                  shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &superior_edge_triangle_map,
                                                                          const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ec) {
        auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto& e1:*superior_edge_set){
            if(evicted_set->count(e1)){
                continue;
            }
            Ec->insert(e1);
            for(const auto& e2:*superior_edge_set){
                if(e1 == e2){
                    continue;
                }

                for (const auto &e3: *superior_edge_triangle_map->at(e1)) {
                    if (superior_edge_triangle_map->at(e2)->count(e3)) {
                        evicted_set->insert(e2);
                        break;
                    }
                }
            }
        }
    }


    uint32_t jes_truss_maintenance::compute_truss_support(const shared_ptr<abstract_graph> &G,
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


    void jes_truss_maintenance::compute_insert_superior_edge_set(const shared_ptr<abstract_graph>&G,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& superior_edge_set,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>,
                                                                         shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>& superior_edge_triangle_map) {
        auto result_map = make_shared<unordered_map<shared_ptr<abstract_edge>,uint32_t>>();
        for (const auto &e1:*edge_set) {
            result_map->insert({e1, 0});
            superior_edge_triangle_map->insert({e1, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            superior_edge_triangle_map->at(e1)->insert(e1);
        }

        for (const auto &e1:*edge_set) {
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
                superior_edge_triangle_map->at(e1)->insert(e2);
                superior_edge_triangle_map->at(e1)->insert(e3);
            }
        }
        for(const auto& [e,superior_support]:*result_map)
        {
            if(superior_support == 0){
                edge_set->erase(e);
            }

            if(superior_support >= 1)
            {
                superior_edge_set->insert(e);
            }
        }
    }

    void jes_truss_maintenance::compute_delete_superior_edge_set(const shared_ptr<abstract_graph>&G,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                                 const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& superior_edge_set,
                                                                 const shared_ptr<unordered_map<shared_ptr<abstract_edge>,
                                                                         shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>& superior_edge_triangle_map) {
        auto result_map = make_shared<unordered_map<shared_ptr<abstract_edge>,uint32_t>>();
        for (const auto &e1:*edge_set) {
            result_map->insert({e1, 0});
            superior_edge_triangle_map->insert({e1, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            superior_edge_triangle_map->at(e1)->insert(e1);
        }
        for (const auto &e1:*edge_set) {
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
                superior_edge_triangle_map->at(e1)->insert(e2);
                superior_edge_triangle_map->at(e1)->insert(e3);
            }
        }

        for(const auto& [e,superior_support]:*result_map)
        {
            if(superior_support == 0){
                edge_set->erase(e);
                edge_truss_map->erase(e);
                G->remove_edge(e);
            }

            if(superior_support >= 1)
            {
                superior_edge_set->insert(e);
            }
        }
    }

    void jes_truss_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &EI,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                       const shared_ptr<thread_pool>& pool) {
        for(const auto &e:*EI){
            edge_truss_map->insert({e, 2});
        }

        while (!EI->empty()) {
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            compute_insert_edge_set(G, EI, edge_truss_map, Ec, EU);

            auto candidate_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            for (const auto &[k,EUk]:*EU) {
                candidate_edge_map->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }

            /**
              * @brief parallel part
              */
            for (const auto &k_pair:*EU) {
                pool->submit_task([=] {
                    auto k = k_pair.first;
                    auto EUk = k_pair.second;
                    auto candidate_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();
                    k_joint_insert(G, EUk, edge_truss_map,candidate_edge_map->at(k), candidate_edge_support_map, k);
                });
            }
            pool->barrier();

            for(const auto& k_pair:*candidate_edge_map)
            {
                for(const auto& e:*k_pair.second)
                {
                    edge_truss_map->at(e) = edge_truss_map->at(e) + 1;
                }
            }

            for (const auto &e:*Ec) {
                auto truss_number = compute_truss_number(G, edge_truss_map, e);
                edge_truss_map->at(e) = truss_number;
            }
        }
    }

    void jes_truss_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                       const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &ED,
                                       const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                       const shared_ptr<thread_pool>& pool) {
        while (!ED->empty()) {
            auto Ec = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            auto EU = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            compute_delete_edge_set(G, ED, edge_truss_map, Ec, EU);

            auto candidate_edge_map = make_shared<map<uint32_t,shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
            for (const auto &[k,EUk]:*EU) {
                candidate_edge_map->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            }

            for(const auto &e:*Ec){
                G->remove_edge(e);
                ED->erase(e);
                edge_truss_map->erase(e);
            }

            /**
             * @brief parallel part
             */
            for (const auto &k_pair:*EU) {
                pool->submit_task([=]() {
                    auto k = k_pair.first;
                    auto Ek = k_pair.second;
                    k_joint_delete(G, Ek, edge_truss_map, candidate_edge_map->at(k), k);
                });
            }
            pool->barrier();

            for(const auto&[k,Ek]:*candidate_edge_map)
            {
                for (const auto &e:*Ek) {
                    --edge_truss_map->at(e);
                }
            }
        }
    }


    void jes_truss_maintenance::k_joint_delete(const shared_ptr<abstract_graph> &G,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &Ek,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                               uint32_t k) {
        auto cd = make_shared<unordered_map<shared_ptr<abstract_edge>, uint32_t>>();

        auto evicted = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for(const auto &e:*Ek){
            cd->insert({e, compute_truss_support(G, edge_truss_map, e, k)});
            if(cd->at(e) < k - 2){
                evicted->insert(e);
            }
        }

        while(!evicted->empty()){
            auto e1 = *evicted->begin();
            evicted->erase(e1);

            candidate_edge_set->insert(e1);

            auto u = e1->get_source_vertex_id();
            auto v = e1->get_destination_vertex_id();

            if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                swap(u, v);
            }

            for(const auto&[w,e2]:*G->get_vertex(u)->get_edge_map()){
                if(w == v || edge_truss_map->at(e2) < k || candidate_edge_set->count(e2)){
                    continue;
                }
                auto e3 = G->get_edge(v,w);
                if(!e3 || edge_truss_map->at(e3) < k || candidate_edge_set->count(e3)){
                    continue;
                }

                if(edge_truss_map->at(e2) == k){
                    if(!cd->count(e2)){
                        cd->insert({e2, compute_truss_support(G, edge_truss_map, e2,k)});
                    }
                    --cd->at(e2);
                    if(cd->at(e2) < k - 2){
                        evicted->insert(e2);
                    }
                }

                if(edge_truss_map->at(e3) == k){
                    if(!cd->count(e3)){
                        cd->insert({e3, compute_truss_support(G, edge_truss_map, e3, k)});
                    }
                    --cd->at(e3);
                    if(cd->at(e3) < k - 2){
                        evicted->insert(e3);
                    }
                }
            }
        }
    }

    void jes_truss_maintenance::k_joint_insert(const shared_ptr<abstract_graph> &G,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &current_edge_set,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                               const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &candidate_edge_set,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                               uint32_t k) {
        auto evicted_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        auto visited_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set){
                candidate_edge_support_map->insert({e1, 0});

                auto e1_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u, v);
                }

                for(const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()){
                    if(w == v || evicted_set->count(e2)|| edge_truss_map->at(e2) < k){
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if(!e3 || evicted_set->count(e3)|| edge_truss_map->at(e3) < k){
                        continue;
                    }

                    ++candidate_edge_support_map->at(e1);

                    if(edge_truss_map->at(e2) == k && !visited_set->count(e2) && !current_edge_set->count(e2)){
                        e1_set->insert(e2);
                    }

                    if(edge_truss_map->at(e3) == k && !visited_set->count(e3) && !current_edge_set->count(e3)){
                        e1_set->insert(e3);
                    }
                }

                if(candidate_edge_support_map->at(e1) > k - 2){
                    candidate_edge_set->insert(e1);
                    next_edge_set->merge(*e1_set);
                }else
                {
                    remove_unsatisfied_edges(G, edge_truss_map, e1, candidate_edge_set, candidate_edge_support_map, evicted_set, k);
                }
            }
            visited_set->merge(*current_edge_set);
            swap(*current_edge_set, *next_edge_set);
        }
    }

    void jes_truss_maintenance::remove_unsatisfied_edges(const shared_ptr<abstract_graph> &G,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &edge_truss_map,
                                                         const shared_ptr<abstract_edge>& e,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& candidate_edge_set,
                                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>> &candidate_edge_support_map,
                                                         const shared_ptr<unordered_set<shared_ptr<abstract_edge>>>& evicted_set,
                                                         uint32_t k)
    {
        auto current_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        current_edge_set->insert(e);

        while(!current_edge_set->empty()){
            auto next_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
            for(const auto &e1:*current_edge_set){
                candidate_edge_set->erase(e1);
                evicted_set->insert(e1);

                auto u = e1->get_source_vertex_id();
                auto v = e1->get_destination_vertex_id();

                if(G->get_vertex(u)->get_degree() > G->get_vertex(v)->get_degree()){
                    swap(u, v);
                }

                for(const auto &[w,e2]:*G->get_vertex(u)->get_edge_map()){
                    if(w == v || edge_truss_map->at(e2) < k || evicted_set->count(e2)){
                        continue;
                    }

                    auto e3 = G->get_edge(v, w);
                    if(!e3 || edge_truss_map->at(e3) < k || evicted_set->count(e3)){
                        continue;
                    }

                    if(candidate_edge_set->count(e2) && candidate_edge_support_map->at(e2) > k -2){
                        --candidate_edge_support_map->at(e2);
                        if(candidate_edge_support_map->at(e2) <= k - 2 && !current_edge_set->count(e2)){
                            next_edge_set->insert(e2);
                        }
                    }

                    if(candidate_edge_set->count(e3)  && candidate_edge_support_map->at(e3) > k -2){
                        --candidate_edge_support_map->at(e3);
                        if(candidate_edge_support_map->at(e3) <= k - 2 && !current_edge_set->count(e3)){
                            next_edge_set->insert(e3);
                        }
                    }
                }
            }
            swap(*current_edge_set, *next_edge_set);
        }
    }
}