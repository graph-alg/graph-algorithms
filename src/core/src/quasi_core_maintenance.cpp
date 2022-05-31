
#include "core/quasi_core_maintenance.h"

namespace scnu {
    /**
    * @details find a candidate graph by using quasi-k-core
    * @param G
    * @param quasi_core_edge_map
    * @param vertex_core_map
    * @param k
    * @return
    */
    void quasi_core_maintenance::candidate_graph_finding(const shared_ptr<abstract_graph> &G,
                                                        const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &quasi_core_edge_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& CD,
                                                        const shared_ptr<unordered_set<uint32_t>>& Vc,
                                                        const shared_ptr<unordered_map<uint32_t, uint32_t>>& Vc_degree_map,
                                                        uint32_t k) {
        /**
         * @brief a seed set
         */
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        auto visited_set = make_shared<unordered_set<uint32_t>>();
        for (const auto &k_pair:*quasi_core_edge_map) {
            if (k_pair.first >= k) {
                for (const auto &e:*k_pair.second) {
                    auto u = e->get_source_vertex_id();
                    auto v = e->get_destination_vertex_id();

                    /**
                     * @brief ignore invalid vertices
                     */
                    if (vertex_core_map->at(u) >= k-1 && vertex_core_map->at(v) >= k-1) {

                        if (vertex_core_map->at(u) == k - 1) {
                            if (CD->at(u) == 0) {
                                auto degree = get_core_degree(G, vertex_core_map, u, k - 1);
                                CD->at(u) = degree;
                            }
                            if (CD->at(u) >= k) {
                                vertex_set->insert(u);
                            }
                        }

                        if (vertex_core_map->at(v) == k - 1) {
                            if (CD->at(v) == 0) {
                                auto degree = get_core_degree(G, vertex_core_map, v, k - 1);
                                CD->at(v) =  degree;
                            }
                            if (CD->at(v) >= k) {
                                vertex_set->insert(v);
                            }
                        }
                    }
                }
            }
        }

        auto evicted_set = make_shared<unordered_set<uint32_t>>();

        while (!vertex_set->empty()) {
            auto u = *vertex_set->begin();
            vertex_set->erase(u);

            if(visited_set->count(u)){
                continue;
            }
            visited_set->insert(u);

            auto u_degree = 0;
            auto u_set = make_shared<unordered_set<uint32_t>>();
            for (const auto &[v, e]: *G->get_vertex(u)->get_edge_map()) {
                if(vertex_core_map->at(v) >= k-1 && !evicted_set->count(v)){
                    if (Vc->count(v) || vertex_core_map->at(v) >= k) {
                        ++u_degree;
                        continue;
                    }

                    if(vertex_core_map->at(v) == k-1){
                        if(CD->at(v) == 0){
                            auto degree = get_core_degree(G, vertex_core_map, v, k - 1);
                            CD->at(v) = degree;
                        }

                        if(CD->at(v) >= k){
                            ++u_degree;
                            u_set->insert(v);
                        }
                    }
                }
            }
            if (u_degree >= k) {
                Vc->insert(u);
                Vc_degree_map->at(u) = u_degree;
                vertex_set->merge(*u_set);
            } else {
                remove_unsatisfied_vertices(G,vertex_core_map, Vc,u, k, Vc_degree_map, evicted_set);
            }
        }

    }

    /**
     * @details find a bottom-up candidate graph by using quasi-k-core
     * @param G
     * @param edge_quasi_core_map
     * @param vertex_core_map
     * @param CD
     * @param k
     * @param thread_number
     * @return
     */
    void quasi_core_maintenance::candidate_graph_finding(const shared_ptr<abstract_graph> &G,
                                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                         const shared_ptr<unordered_map<uint32_t, uint32_t>> & CD,
                                                         const shared_ptr<unordered_set<uint32_t>> &V,
                                                         const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                         const shared_ptr<unordered_map<uint32_t ,uint32_t>> &Vc_degree_map,
                                                         const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& current_task_vector,
                                                         const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& evicted_task_vector,
                                                         uint32_t k,
                                                         const shared_ptr<thread_pool> &pool) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        {
         /**
           * @brief a seed set
           */
            for (auto iter = V->begin();iter!=V->end();) {
                auto v = *iter;
                ++iter;
                if (vertex_core_map->at(v) < k - 1) {
                    V->erase(v);
                } else if (vertex_core_map->at(v) == k - 1) {
                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        CD->at(v) = get_core_degree(G, vertex_core_map, v, k - 1);
                        if (CD->at(v) >= k) {
                            current_task_vector->at(thread_id)->insert(v);
                        }
                    });
                }
            }
            pool->barrier();
            merger_set(current_task_vector, vertex_set, pool);
        }

        auto visited_set = make_shared<unordered_set<uint32_t>>();
        {
            while(!vertex_set->empty()){
                for (const auto&u:*vertex_set) {
                    if(visited_set->count(u)){
                        continue;
                    }
                    visited_set->insert(u);
                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());

                        auto u_degree = 0;
                        auto u_set = make_shared<unordered_set<uint32_t>>();
                        for (const auto &p: *G->get_vertex(u)->get_edge_map()) {
                            auto v = p.first;

                            if (vertex_core_map->at(v) >= k - 1) {
                                if(Vc->count(v) || vertex_core_map->at(v) >= k){
                                    ++u_degree;
                                    continue;
                                }

                                if(vertex_core_map->at(v) == k-1){
                                    vertex_mutex_map->at(v)->lock();
                                    if(CD->at(v) == 0){
                                        CD->at(v) = get_core_degree(G, vertex_core_map, v, k - 1);
                                    }
                                    vertex_mutex_map->at(v)->unlock();


                                    if (CD->at(v) >= k) {
                                         ++u_degree;
                                         u_set->insert(v);
                                    }
                                }
                            }
                        }

                        if(u_degree >= k){
                            Vc_degree_map->at(u) = u_degree;
                            current_task_vector->at(thread_id)->merge(*u_set);
                        }else
                        {
                            evicted_task_vector->at(thread_id)->insert(u);
                        }
                   });
                }
                pool->barrier();
                Vc->merge(*vertex_set);
                vertex_set->clear();
                merger_set(current_task_vector, vertex_set, pool);
            }
            remove_unsatisfied_vertices(G, vertex_mutex_map, vertex_core_map, Vc, Vc_degree_map, evicted_task_vector, k, pool);
        }
    }

    /**
     * @details insert maintenance based on graph view
     * @param previous_graph
     * @param insertion_edge_set
     * @param vertex_core_map
     * @param k_max
     */
    void quasi_core_maintenance::insert(const shared_ptr<abstract_graph> &previous_graph,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &insertion_edge_set,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                        const shared_ptr<uint32_t>& k_max) {
        /**
         * @brief insert edges
         */
        for (const auto &e: *insertion_edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();
            previous_graph->insert_edge(e);
            if (!vertex_core_map->count(u)) {
                vertex_core_map->insert({u, 1});
            }
            if (!vertex_core_map->count(v)) {
                vertex_core_map->insert({v, 1});
            }
        }

        /**
         * @brief quasi core decomposition for edges
         */
        auto quasi_core_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        auto quasi_max_k = quasi_core_edge_decomposition(previous_graph, insertion_edge_set, quasi_core_edge_map);

        auto Vc = make_shared<unordered_set<uint32_t>>();
        auto CD = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto Vc_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();

        for(const auto&[u, core_number]:*vertex_core_map)
        {
            CD->insert({u,0});
            Vc_degree_map->insert({u,0});
        }

        uint32_t k = 2;

        while (k <= quasi_max_k) {
            if (k <= *k_max) {
                candidate_graph_finding(previous_graph, quasi_core_edge_map, vertex_core_map, CD, Vc, Vc_degree_map, k);
                for(const auto &u:*Vc){
                    vertex_core_map->at(u) = k;
                    Vc_degree_map->at(u) = 0;
                    CD->at(u) = 0;
                }
                Vc->clear();
            } else {
                candidate_graph_finding(previous_graph, quasi_core_edge_map, vertex_core_map, CD, Vc, Vc_degree_map, k);
                auto new_k_max = partial_core_decomposition(previous_graph, Vc, vertex_core_map, Vc_degree_map, k);
                if (new_k_max > 0) {
                    *k_max = new_k_max;
                }
                break;
            }
            ++k;
        }
    }

    /**
     * @details insert maintenance
     * @param G
     * @param edge_set
     * @param vertex_core_map
     * @param k_max
     * @param thread_number
     */
    void
    quasi_core_maintenance::insert(const shared_ptr<abstract_graph> &G,
                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                   const shared_ptr<uint32_t> &k_max,
                                   uint32_t thread_number) {
        auto V = make_shared<unordered_set<uint32_t>>();
        for (const auto &e: *edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            if (!vertex_core_map->count(u)) {
                vertex_core_map->insert({u, 1});
            }

            if (!vertex_core_map->count(v)) {
                vertex_core_map->insert({v, 1});
            }

            G->insert_edge(e);

            V->insert(u);
            V->insert(v);
        }


        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>();
        auto evicted_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>();
        for (uint32_t i = 0; i < thread_number; ++i) {
            current_task_vector->emplace_back(make_shared<unordered_set<uint32_t>>());
            evicted_task_vector->emplace_back(make_shared<unordered_set<uint32_t>>());
        }

        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
        auto Vc = make_shared<unordered_set<uint32_t>>();
        auto Vc_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto CD = make_shared<unordered_map<uint32_t,uint32_t>>();
        for(const auto&[u, core_number]:*vertex_core_map)
        {
            vertex_mutex_map->insert({u, make_shared<mutex>()});
            CD->insert({u,0});
            Vc_degree_map->insert({u,0});
        }

        auto pool = make_shared<thread_pool>(thread_number);

        uint32_t k = 2;
        while (!V->empty()) {
            if(k <= *k_max){
                candidate_graph_finding(G, vertex_mutex_map, vertex_core_map, CD, V, Vc, Vc_degree_map,
                                        current_task_vector, evicted_task_vector, k, pool);

                for (const auto &u: *Vc) {
                    vertex_core_map->at(u) = k;
                    Vc_degree_map->at(u) = 0;
                    CD->at(u) = 0;
                }
                Vc->clear();
                ++k;
            }
            else{
                candidate_graph_finding(G, vertex_mutex_map, vertex_core_map, CD, V, Vc, Vc_degree_map,
                                        current_task_vector, evicted_task_vector, k, pool);
                auto new_k_max = partial_core_decomposition(G, vertex_mutex_map, vertex_core_map, Vc, Vc_degree_map, current_task_vector, k, pool);
                if(new_k_max > 0){
                    *k_max = new_k_max;
                }
                break;
            }

        }
    }

    void quasi_core_maintenance::remove(const shared_ptr<abstract_graph> &previous_graph,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &removal_edge_set,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                        const shared_ptr<uint32_t>& k_max) {
        /**
         * @brief decompose removed edges to a map of edges and quasi core numbers
         */
        auto quasi_core_edge_map = make_shared<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>();
        auto quasi_max_k = quasi_core_edge_decomposition(previous_graph, removal_edge_set, quasi_core_edge_map);

        /**
         * @remarks when quasi max k is larger than max k, then it is invalid
         */
        if (*k_max < quasi_max_k) {
            quasi_max_k = *k_max;
        }

        auto current_vertex_core_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto CD = make_shared<unordered_map<uint32_t,uint32_t>>();
        for(const auto &[u,core_number]:*vertex_core_map){
            current_vertex_core_map->insert({u, core_number});
        }

        for(const auto&e:*removal_edge_set){
            previous_graph->remove_edge(e);
        }


        uint32_t k = 1;
        auto previous_removed_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto current_removed_vertex_set = make_shared<unordered_set<uint32_t>>();
        while (k <= quasi_max_k) {
            update_single_core(previous_graph, vertex_core_map, quasi_core_edge_map, CD, previous_removed_vertex_set,current_removed_vertex_set, k);
            quasi_core_edge_map->erase(k);
            CD->clear();
            if(k == 1){
                for(const auto&u:*current_removed_vertex_set){
                    current_vertex_core_map->erase(u);
                }
            }else{
                for (const auto &u: *current_removed_vertex_set) {
                    if (!previous_removed_vertex_set->count(u)) {
                        current_vertex_core_map->at(u) = k - 1;
                    }
                }
            }
            previous_removed_vertex_set->clear();
            swap(*previous_removed_vertex_set,*current_removed_vertex_set);
            ++k;
        }

        vertex_core_map->clear();
        swap(*vertex_core_map,*current_vertex_core_map);

        uint32_t max_core_number = 0;
        for(const auto&[v,core_number]:*vertex_core_map){
            if(core_number > max_core_number){
                max_core_number = core_number;
            }
        }
        *k_max = max_core_number;
    }

//    void quasi_core_maintenance::remove(const shared_ptr<abstract_graph> &G,
//                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &removal_edge_set,
//                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
//                                        const shared_ptr<uint32_t> &k_max,
//                                        uint32_t thread_number) {
//        auto current_vertex_core_map = make_shared<unordered_map<uint32_t,uint32_t>>();
//        auto CD = make_shared<unordered_map<uint32_t, uint32_t>>();
//        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();
//
//        for (const auto &[u,core_number]:*vertex_core_map) {
//            CD->insert({u, 0});
//            vertex_mutex_map->insert({u, make_shared<mutex>()});
//            current_vertex_core_map->insert({u, core_number});
//        }
//
//        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
//        auto V = make_shared<unordered_set<uint32_t>>();
//
//        uint32_t k = 0;
//        for (const auto &e: *removal_edge_set) {
//            auto u = e->get_source_vertex_id();
//            auto v = e->get_destination_vertex_id();
//
//            auto min_k = std::min(vertex_core_map->at(u), vertex_core_map->at(v));
//            if(min_k > k){
//                k = min_k;
//            }
//            G->remove_edge(e,isolated_vertex_set);
//            V->insert(u);
//            V->insert(v);
//        }
//
//        for (const auto &u: *isolated_vertex_set) {
//            V->erase(u);
//            current_vertex_core_map->erase(u);
//        }
//
//        bool update_flag = false;
//        if(k == *k_max){
//            update_flag = true;
//        }
//
//        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>();
//        for(uint32_t i = 0;i<thread_number;++i){
//            current_task_vector->emplace_back(make_shared<unordered_set<uint32_t>>());
//        }
//
//        auto current_removed_vertex_set = make_shared<unordered_set<uint32_t>>();
//        auto pool = make_shared<thread_pool>(thread_number);
//        while (k > 0) {
//            update_single_core(G, vertex_mutex_map, vertex_core_map, current_vertex_core_map,
//                               current_removed_vertex_set, CD, V, k, current_task_vector, pool);
//            for(const auto &u:*current_removed_vertex_set){
//                current_vertex_core_map->at(u) = k - 1;
//                CD->at(u) = 0;
//            }
//            current_removed_vertex_set->clear();
//            --k;
//        }
//
//        vertex_core_map->clear();
//        swap(*vertex_core_map,*current_vertex_core_map);
//
//        if(update_flag){
//            uint32_t max_core_number = 0;
//            for(const auto&[v,core_number]:*vertex_core_map){
//                if(core_number > max_core_number){
//                    max_core_number = core_number;
//                }
//            }
//            *k_max = max_core_number;
//        }
//
//    }

    void quasi_core_maintenance::remove(const shared_ptr<abstract_graph> &G,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &removal_edge_set,
                                        const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                        const shared_ptr<uint32_t> &k_max,
                                        uint32_t thread_number) {
        auto CD = make_shared<unordered_map<uint32_t, uint32_t>>();
        auto vertex_mutex_map = make_shared<unordered_map<uint32_t, shared_ptr<mutex>>>();

        for (const auto &[u,core_number]:*vertex_core_map) {
            CD->insert({u, 0});
            vertex_mutex_map->insert({u, make_shared<mutex>()});
        }

        auto isolated_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto V = make_shared<unordered_set<uint32_t>>();

        uint32_t k = 0;
        for (const auto &e: *removal_edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            auto min_k = std::min(vertex_core_map->at(u), vertex_core_map->at(v));
            if(min_k > k){
                k = min_k;
            }
            G->remove_edge(e,isolated_vertex_set);
            V->insert(u);
            V->insert(v);
        }

        for (const auto &u: *isolated_vertex_set) {
            V->erase(u);
            vertex_core_map->erase(u);
        }

        bool update_flag = false;
        if(k == *k_max){
            update_flag = true;
        }

        auto current_task_vector = make_shared<vector<shared_ptr<unordered_set<uint32_t>>>>();
        for(uint32_t i = 0;i<thread_number;++i){
            current_task_vector->emplace_back(make_shared<unordered_set<uint32_t>>());
        }

        auto current_removed_vertex_set = make_shared<unordered_set<uint32_t>>();
        auto pool = make_shared<thread_pool>(thread_number);
        while (k > 0) {
            update_single_core(G, vertex_mutex_map, vertex_core_map,
                               current_removed_vertex_set, CD, V, k, current_task_vector, pool);
            for(const auto &u:*current_removed_vertex_set){
                vertex_core_map->at(u) = k - 1;
                CD->at(u) = 0;
            }
            current_removed_vertex_set->clear();
            --k;
        }

        if(update_flag){
            uint32_t max_core_number = 0;
            for(const auto&[v,core_number]:*vertex_core_map){
                if(core_number > max_core_number){
                    max_core_number = core_number;
                }
            }
            *k_max = max_core_number;
        }

    }


    uint32_t quasi_core_maintenance::get_core_degree(const shared_ptr<abstract_graph> &G,
                                                     const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                     uint32_t u,
                                                     uint32_t k) {
        uint32_t core_degree = 0;
        auto u_vertex = G->get_vertex(u);
        if(!u_vertex){
            return core_degree;
        }
        for(const auto &[v,e]:*u_vertex->get_edge_map()){
            if(vertex_core_map->at(v) >= k){
                ++core_degree;
            }
        }
        return core_degree;
    }

    void quasi_core_maintenance::merger_set(const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &input_set_vector,
                                            const shared_ptr<unordered_set<uint32_t>> &output_set,
                                            const shared_ptr<thread_pool> &pool){
        for(auto i = 1; i < input_set_vector->size();i *= 2){
            for(uint32_t j = 0; j + i < input_set_vector->size();j += i * 2){
                pool->submit_task([=]{
                    auto set1 = input_set_vector->at(j);
                    auto set2 = input_set_vector->at(j + i);
                    set1->merge(*set2);
                    set2->clear();
               });
            }
            pool->barrier();
        }
        swap(*output_set,*input_set_vector->at(0));
    }



    /**
     * @details decompose the remaining graph to a map of vertices and core numbers
     * @param G
     * @param Vc
     * @param vertex_core_map
     * @param CD
     * @param k
     * @return
     */
    uint32_t quasi_core_maintenance::partial_core_decomposition(
            const shared_ptr<abstract_graph> &G,
            const shared_ptr<unordered_set<uint32_t>> &Vc,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
            const shared_ptr<unordered_map<uint32_t, uint32_t>> &CD,
            uint32_t k) {
        uint32_t k_max = 0;
        auto evicted_set = make_shared<unordered_set<uint32_t>>();
        while (!Vc->empty()) {
            k_max = k;

            for(const auto &v:*Vc){
                if(CD->at(v) <= k){
                    evicted_set->insert(v);
                }
            }

            while(!evicted_set->empty()){
                auto u = *evicted_set->begin();
                evicted_set->erase(u);
                Vc->erase(u);
                vertex_core_map->at(u) = k;

                for(const auto &[v,e]:*G->get_vertex(u)->get_edge_map()){
                    if(!Vc->count(v)){
                        continue;
                    }
                    --CD->at(v);
                    if(CD->at(v) <= k){
                        evicted_set->insert(v);
                    }
                }
            }
            ++k;
        }
        return k_max;
    }

    /**
     * @details decompose the remaining graph
     * @param G
     * @param Vc
     * @param vertex_core_map
     * @param k
     * @param thread_number
     * @return
     */
    uint32_t quasi_core_maintenance::partial_core_decomposition(const shared_ptr<abstract_graph> &G,
                                                                const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vector_mutex_vector,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                                const shared_ptr<unordered_set<uint32_t>> &F_k,
                                                                const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_degree_map,
                                                                const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>> &current_task_vector,
                                                                uint32_t k,
                                                                const shared_ptr<thread_pool>& pool) {
        uint32_t k_max = 0;
        while (!F_k->empty()) {
            k_max = k;
            pool->submit_task([=]{
                for (const auto &u:*F_k) {
                    if (vertex_degree_map->at(u) < k + 1) {
                        current_task_vector->at(0)->insert(u);
                    }
                }
            },0 % pool->get_thread_number());
            for (const auto &u:*F_k) {
                vertex_core_map->at(u) = k;
            }
            pool->barrier();


            remove_unsatisfied_vertices(G, vector_mutex_vector, vertex_core_map, F_k, vertex_degree_map,
                                        current_task_vector, k + 1, pool);
            ++k;
        }
        return k_max;
    }

    uint32_t quasi_core_maintenance::quasi_core_edge_decomposition(const shared_ptr<abstract_graph> &entire_graph,
                                                                   const shared_ptr<unordered_set<shared_ptr<abstract_edge>>> &edge_set,
                                                                   const shared_ptr<unordered_map<uint32_t,
                                                                           shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>>
                                                                   &quasi_core_edge_map) {
        /**
         * @details get the global degree (entire graph) of vertices of inserted edges
         */
        auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        /**
         * @details get the neighbor set of vertices of inserted edges
         * @remarks the size of neighbor set is not equal to degree
         */
        auto sub_graph = make_shared<abstract_graph>();

        auto evicted_edge_set = make_shared<unordered_set<shared_ptr<abstract_edge>>>();
        for (const auto &e:*edge_set) {
            auto u = e->get_source_vertex_id();
            auto v = e->get_destination_vertex_id();

            evicted_edge_set->insert(e);

            sub_graph->insert_edge(e);

            if (!vertex_degree_map->count(u)) {
                auto u_vertex = entire_graph->get_vertex(u);
                auto degree = u_vertex->get_degree();
                vertex_degree_map->insert({u, degree});
            }

            if (!vertex_degree_map->count(v)) {
                auto v_vertex = entire_graph->get_vertex(v);
                auto degree = v_vertex->get_degree();
                vertex_degree_map->insert({v, degree});
            }
        }

        uint32_t k_max = 1;
        uint32_t k = 1;

        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        while (!evicted_edge_set->empty()) {
            for (const auto &[u, sub_u_vertex]:*sub_graph->get_vertex_map()) {
                if (vertex_degree_map->at(u) <= k) {
                    vertex_set->insert(u);
                }
            }

            quasi_core_edge_map->insert({k, make_shared<unordered_set<shared_ptr<abstract_edge>>>()});
            k_max = k;


            while (!vertex_set->empty()) {
                auto u = *vertex_set->begin();
                vertex_set->erase(u);

                /**
                 * @remarks in insert graph rather than global graph
                 */
                for (const auto &[v,e]:*sub_graph->get_vertex(u)->get_edge_map()) {
                    sub_graph->get_vertex(v)->remove_edge(u);

                    /**
                     * @brief set the quasi core number for edge
                     */
                    quasi_core_edge_map->at(k)->insert(e);
                    evicted_edge_set->erase(e);

                    --vertex_degree_map->at(v);
                    if (vertex_degree_map->at(v) <= k) {
                        vertex_set->insert(v);
                    }
                }

                sub_graph->remove_vertex(u);
            }
            ++k;
        }
        return k_max;
    }

    /**
    * @details recursively remove unsatisfied vertices
    * @param G
    * @param Vc
    * @param w
    * @param k
    * @param Vc_degree_map
    * @param evicted_set
    */
    void quasi_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                             const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                             uint32_t w,
                                                             uint32_t k,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &Vc_degree_map,
                                                             const shared_ptr<unordered_set<uint32_t>>& evicted_set) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        vertex_set->insert(w);

        while (!vertex_set->empty()) {
            auto u = *vertex_set->begin();
            vertex_set->erase(u);
            evicted_set->insert(u);

            auto u_vertex = G->get_vertex(u);
            for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                if (vertex_core_map->at(v) == k-1 && Vc_degree_map->at(v) >= k) {
                    --Vc_degree_map->at(v);
                    if (Vc_degree_map->at(v) < k) {
                        vertex_set->insert(v);
                    }
                }
            }

            Vc->erase(u);
            Vc_degree_map->at(u) = 0;
        }
    }

    void quasi_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                             const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                             const shared_ptr<unordered_set<uint32_t>> &vertex_set,
                                                             uint32_t k,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &Vc_degree_map) {
        while (!vertex_set->empty()) {
            auto u = *vertex_set->begin();
            vertex_set->erase(u);

            auto u_vertex = G->get_vertex(u);
            for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                if (vertex_core_map->at(v) == k-1 && Vc_degree_map->at(v) >= k) {
                    --Vc_degree_map->at(v);
                    if (Vc_degree_map->at(v) < k) {
                        vertex_set->insert(v);
                    }
                }
            }

            Vc->erase(u);
            Vc_degree_map->at(u) = 0;
        }
    }



    void quasi_core_maintenance::remove_unsatisfied_vertices(const shared_ptr<abstract_graph> &G,
                                                             const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_map,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                             const shared_ptr<unordered_set<uint32_t>> &Vc,
                                                             const shared_ptr<unordered_map<uint32_t, uint32_t>> &Vc_degree_map,
                                                             const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& current_task_vector,
                                                             uint32_t k,
                                                             const shared_ptr<thread_pool>& pool)
    {
        auto current_task_set = make_shared<unordered_set<uint32_t>>();
        merger_set(current_task_vector, current_task_set, pool);

        while (!current_task_set->empty()){
            for(const auto& u:*current_task_set) {
                Vc->erase(u);
                Vc_degree_map->at(u) = 0;
                pool->submit_task([=] {
                    auto thread_id = pool->get_thread_id(std::this_thread::get_id());

                    auto u_vertex = G->get_vertex(u);
                    for (const auto &p:*u_vertex->get_edge_map()) {
                        auto v = p.first;
                        if(vertex_core_map->at(v) == k-1 && Vc_degree_map->at(v) >= k){

                            vertex_mutex_map->at(v)->lock();
                            --Vc_degree_map->at(v);
                            vertex_mutex_map->at(v)->unlock();

                            if (Vc_degree_map->at(v) < k) {
                                current_task_vector->at(thread_id)->insert(v);
                            }
                        }
                    }
                });
            }
            current_task_set->clear();
            pool->barrier();
            merger_set(current_task_vector, current_task_set, pool);
        }
    }


    /**
     * @details update the k-core by recursively deleting removed edges
     * @param G
     * @param quasi_core_edge_map
     * @param k

     */
    void quasi_core_maintenance::update_single_core(const shared_ptr<abstract_graph> &G,
                                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                                                   const shared_ptr<unordered_map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_edge>>>>> &quasi_core_edge_map,
                                                                                   const shared_ptr<unordered_map<uint32_t, uint32_t>>& CD,
                                                                                   const shared_ptr<unordered_set<uint32_t>>& previous_removed_vertex_set,
                                                                                   const shared_ptr<unordered_set<uint32_t>>& current_removed_vertex_set,
                                                                                   uint32_t k) {
        auto vertex_set = make_shared<unordered_set<uint32_t>>();
        for (const auto &k_pair:*quasi_core_edge_map) {
            if (k_pair.first >= k) {
                for (const auto &e:*k_pair.second) {
                    auto u = e->get_source_vertex_id();
                    auto v = e->get_destination_vertex_id();

                    if (vertex_core_map->at(u) >= k &&  vertex_core_map->at(v) >= k) {
                        if(previous_removed_vertex_set->count(u)){
                            vertex_set->insert(u);
                        }
                        else
                        {
                            if(!CD->count(u)){
                                auto degree = get_core_degree(G, vertex_core_map, u, k);
                                CD->insert({u,degree});
                            }

                            if (CD->at(u) < k) {
                                vertex_set->insert(u);
                            }
                        }

                        if(previous_removed_vertex_set->count(v)){
                            vertex_set->insert(v);
                        }
                        else
                        {
                            if(!CD->count(v)){
                                auto degree = get_core_degree(G, vertex_core_map, v, k);
                                CD->insert({v,degree});
                            }

                            if (CD->at(v) < k) {
                                vertex_set->insert(v);
                            }
                        }
                    }
                }
            }
        }

        while (!vertex_set->empty()) {
            auto u = *vertex_set->begin();
            vertex_set->erase(u);
            current_removed_vertex_set->insert(u);

            auto u_vertex = G->get_vertex(u);
            if(!u_vertex){
                continue;
            }

            for (const auto &[v, e]:*u_vertex->get_edge_map()) {
                if (vertex_core_map->at(v) >= k) {
                    if(previous_removed_vertex_set->count(v)){
                        if(!current_removed_vertex_set->count(v)){
                            vertex_set->insert(v);
                        }
                    }
                    else
                    {
                        if(!CD->count(v)){
                            auto degree = get_core_degree(G, vertex_core_map, v, k);
                            CD->insert({v, degree});
                        }

                        --CD->at(v);
                        if (CD->at(v) < k && !current_removed_vertex_set->count(v)) {
                            vertex_set->insert(v);
                        }
                    }
                }
            }
        }
    }


//    void quasi_core_maintenance::update_single_core(const shared_ptr<abstract_graph> &G,
//                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_vector,
//                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
//                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &current_vertex_core_map,
//                                                    const shared_ptr<unordered_set<uint32_t>>& current_removed_vertex_set,
//                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &CD,
//                                                    const shared_ptr<unordered_set<uint32_t>> &V,
//                                                    uint32_t k,
//                                                    const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& current_task_vector,
//                                                    const shared_ptr<thread_pool> &pool) {
//        auto current_task_set = make_shared<unordered_set<uint32_t>>();
//        {
//            for (auto iter = V->begin();iter!=V->end();) {
//                auto v = *iter;
//                ++iter;
//                if (current_vertex_core_map->at(v) > k){
//                    V->erase(v);
//                }
//                else if (current_vertex_core_map->at(v) == k) {
//                    pool->submit_task([=] {
//                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
//                        CD->at(v) = get_core_degree(G, vertex_core_map, v, k);
//                        if (CD->at(v) < k) {
//                            current_task_vector->at(thread_id)->insert(v);
//                        }
//                    });
//                }
//            }
//            pool->barrier();
//            merger_set(current_task_vector,current_task_set,pool);
//        }
//
//
//
//        {
//            while (!current_task_set->empty()) {
//                for (const auto&u:*current_task_set) {
//                    pool->submit_task([=]{
//                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
//                        for (const auto &p:*G->get_vertex(u)->get_edge_map()) {
//                            auto v = p.first;
//                            if (vertex_core_map->at(v) >= k && !current_removed_vertex_set->count(v) && !current_task_set->count(v)) {
//                                if(current_vertex_core_map->at(v) == k){
//                                    vertex_mutex_vector->at(v)->lock();
//                                    if (CD->at(v) == 0) {
//                                        CD->at(v) = get_core_degree(G, vertex_core_map, v, k);
//                                    }
//                                    --CD->at(v);
//                                    vertex_mutex_vector->at(v)->unlock();
//
//                                    if (CD->at(v) < k) {
//                                        current_task_vector->at(thread_id)->insert(v);
//                                    }
//                                }
//                            }
//                        }
//                    });
//                }
//                pool->barrier();
//                current_removed_vertex_set->merge(*current_task_set);
//                merger_set(current_task_vector, current_task_set, pool);
//            }
//        }
//    }

    void quasi_core_maintenance::update_single_core(const shared_ptr<abstract_graph> &G,
                                                    const shared_ptr<unordered_map<uint32_t, shared_ptr<mutex>>>& vertex_mutex_vector,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_core_map,
                                                    const shared_ptr<unordered_set<uint32_t>>& current_removed_vertex_set,
                                                    const shared_ptr<unordered_map<uint32_t, uint32_t>> &CD,
                                                    const shared_ptr<unordered_set<uint32_t>> &V,
                                                    uint32_t k,
                                                    const shared_ptr<vector<shared_ptr<unordered_set<uint32_t>>>>& current_task_vector,
                                                    const shared_ptr<thread_pool> &pool) {
        auto current_task_set = make_shared<unordered_set<uint32_t>>();
        {
            for (auto iter = V->begin();iter!=V->end();) {
                auto v = *iter;
                ++iter;
                if (vertex_core_map->at(v) > k){
                    V->erase(v);
                }
                else if (vertex_core_map->at(v) == k) {
                    pool->submit_task([=] {
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        CD->at(v) = get_core_degree(G, vertex_core_map, v, k);
                        if (CD->at(v) < k) {
                            current_task_vector->at(thread_id)->insert(v);
                        }
                    });
                }
            }
            pool->barrier();
            merger_set(current_task_vector,current_task_set,pool);
        }



        {
            while (!current_task_set->empty()) {
                for (const auto&u:*current_task_set) {
                    pool->submit_task([=]{
                        auto thread_id = pool->get_thread_id(std::this_thread::get_id());
                        for (const auto &p:*G->get_vertex(u)->get_edge_map()) {
                            auto v = p.first;
                            if (vertex_core_map->at(v) == k && !current_removed_vertex_set->count(v) && !current_task_set->count(v)) {
                                vertex_mutex_vector->at(v)->lock();
                                if (CD->at(v) == 0) {
                                    CD->at(v) = get_core_degree(G, vertex_core_map, v, k);
                                }
                                --CD->at(v);
                                vertex_mutex_vector->at(v)->unlock();

                                if (CD->at(v) < k) {
                                    current_task_vector->at(thread_id)->insert(v);
                                }
                            }
                        }
                    });
                }
                pool->barrier();
                current_removed_vertex_set->merge(*current_task_set);
                merger_set(current_task_vector, current_task_set, pool);
            }
        }
    }
}

