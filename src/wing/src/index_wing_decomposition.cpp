
#include "wing/index_wing_decomposition.h"

namespace scnu{

    /**
     * @details construct sub-bloom-index of epsilon-wing
     * @param G
     * @param edge_set
     * @param edge_support_map
     * @param edge_wing_map
     * @param edge_rank_map
     * @param edge_mutex_map
     * @param bloom_index_mutex
     * @param thread_count
     * @return
     */
    shared_ptr<BE_index> index_wing_decomposition::compressed_index_construction(
            const shared_ptr<abstract_bipartite_graph> &G,
            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,uint32_t>>& edge_support_map,
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
            const shared_ptr<thread_pool>& pool) {

        auto bloom_index = make_shared<BE_index>();
        for(auto iter = edge_set->begin();iter!=edge_set->end();){
            auto &e = *iter;
            ++iter;
            if(edge_wing_map->at(e) >0)
            {
                edge_set->erase(e);
                edge_support_map->erase(e);
            }else
            {
                bloom_index->insert_edge(e,edge_support_map->at(e));
            }
        }
        auto vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        vertex_priority_computation(G,vertex_priority_map, pool);
        left_index_construction(G,vertex_priority_map,edge_support_map,edge_wing_map,edge_mutex_map,bloom_index,pool);
        right_index_construction(G,vertex_priority_map,edge_support_map,edge_wing_map,edge_mutex_map,bloom_index,pool);
        return bloom_index;
    }

    /**
     * @details an optimized decomposition
     * @param G
     * @param edge_wing_map
     * @param thread_number
     */
    void index_wing_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             const shared_ptr<thread_pool>& pool){

        /**
         * @brief initialize
         */
        auto edge_set = G->get_edge_set();
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        /**
         * @brief construct BE_index
         */
        auto bloom_index = index_construction(G, edge_set, edge_support_map, edge_wing_map, edge_mutex_map, pool);

        auto B_mutex_map = make_shared<unordered_map<shared_ptr<priority_obeyed_bloom>,shared_ptr<mutex>>>();
        auto C = make_shared<unordered_map<shared_ptr<priority_obeyed_bloom>,uint32_t>>();
        {
            auto bloom_map = bloom_index->get_bloom_map();
            for(const auto &[p,B]:*bloom_map){
                B_mutex_map->insert({B, shared_ptr<mutex>()});
            }
            auto location_vector = pool->split_task(bloom_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &B = iter->second;
                        B_mutex_map->at(B) =  make_shared<mutex>();
                    }
                });
            }
            for(const auto &[p,B]:*bloom_index->get_bloom_map()){
                C->insert({B, 0});
            }
            pool->barrier();
        }
        pool->barrier();

        auto rank_id = make_shared<uint32_t>(0);
        while (!edge_set->empty()){
            auto MBS = minimal_butterfly_support_finding(bloom_index->get_edge_support_map(),edge_wing_map, pool);
            auto S_cur = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            scan(edge_set, bloom_index->get_edge_support_map(), edge_wing_map, MBS,S_cur, pool);

            auto S_next = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            while (!S_cur->empty()){
                {
                    auto location_vector = pool->split_task(S_cur);
                    for(uint32_t i = 0; i< thread_number; ++i) {
                        pool->submit_task([=] {
                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                auto &e1 = *iter;
                                edge_wing_map->at(e1) = bloom_index->get_support(e1);
                            }
                        });
                    }
                    pool->barrier();
                }

                auto affected_B_set = make_shared<unordered_set<shared_ptr<priority_obeyed_bloom>>>();
                {
                    set_edge_rank(S_cur, edge_rank_map, rank_id);
                    auto location_vector = pool->split_task(S_cur);
                    for(uint32_t i = 0; i< thread_number; ++i){
                        pool->submit_task([=]{
                            auto sub_affected_B_set = make_shared<unordered_set<shared_ptr<priority_obeyed_bloom>>>();

                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &e1 = * iter;
                                //edge_mutex_map->at(e1)->lock();
                                auto B_set = bloom_index->get_bloom_set(e1);
                                for(const auto &B:*B_set){
                                    if(B->get_butterfly_count()>0){

                                        B_mutex_map->at(B)->lock();
                                        auto e2 = B->get_twin(e1);
                                        B_mutex_map->at(B)->unlock();

                                        if(S_cur->count(e2) && edge_rank_map->at(e1) < edge_rank_map->at(e2)){
                                            continue;
                                        }

                                        sub_affected_B_set->insert(B);

                                        B_mutex_map->at(B)->lock();
                                        C->at(B) = C->at(B) + 1;
                                        B_mutex_map->at(B)->unlock();

                                        auto k = B->get_k();

                                        if(edge_wing_map->at(e2) == 0){
                                            edge_mutex_map->at(e2)->lock();
                                            bloom_index->update_support (e2, max(MBS,bloom_index->get_support(e2)-(k-1))) ;
                                            bloom_index->remove_edge(e2,B);
                                            edge_mutex_map->at(e2)->unlock();

                                            B_mutex_map->at(B)->lock();
                                            B->remove_edge(e2);
                                            B_mutex_map->at(B)->unlock();
                                        }
                                    }
                                }
                                //edge_mutex_map->at(e1)->unlock();
                            }

                            global_mutex->lock();
                            affected_B_set->merge(*sub_affected_B_set);
                            global_mutex->unlock();
                        });
                    }
                    pool->barrier();
                }

                {
                    auto location_vector = pool->split_task(affected_B_set);
                    for(uint32_t i = 0; i< thread_number; ++i) {
                        pool->submit_task([=] {
                            auto sub_S_next = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                auto &B = *iter;
                                auto k = B->get_k();
                                B->set_butterfly_count((k - C->at(B)) * ((k - 1) - C->at(B)) / 2);
                                if (B->get_butterfly_count() == 0) {
                                    /**
                                     * @details each thread has different B, so there is no data race
                                     */
                                    C->at(B) = k - 1;
                                }

                                for (const auto&[e1, e2]:*B->get_edge_map()) {
                                    if (!S_cur->count(e1)) {
                                        edge_mutex_map->at(e1)->lock();
                                        bloom_index->update_support(e1,
                                                                    max(MBS, bloom_index->get_support(e1) - C->at(B)));
                                        edge_mutex_map->at(e1)->unlock();

                                        if (bloom_index->get_support(e1) == MBS) {
                                            sub_S_next->insert(e1);
                                        }
                                    }
                                }
                            }

                            global_mutex->lock();
                            S_next->merge(*sub_S_next);
                            global_mutex->unlock();
                        });
                    }
                    pool->barrier();
                }

                {
                    auto location_vector = pool->split_task(affected_B_set);
                    for(uint32_t i = 0; i< thread_number; ++i) {
                        pool->submit_task([=] {
                            auto sub_S_next = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                            auto &sub_begin = *location_vector->at(i);
                            auto &sub_end = *location_vector->at(i + 1);

                            for(auto iter = sub_begin; iter!=sub_end; ++iter) {
                                auto &B = *iter;
                                C->at(B) = 0;
                            }

                        });
                    }
                    for(const auto&e:*S_cur)
                    {
                        edge_set->erase(e);
                        bloom_index->remove_edge(e);
                        edge_rank_map->at(e) = UINT32_MAX;
                    }
                    pool->barrier();
                }
                affected_B_set->clear();
                S_cur->clear();
                swap(*S_cur,*S_next);
            }
        }
    }

    /**
     * @details a top-down decomposition
     * @param G
     * @param edge_wing_map
     * @param tau
     * @param thread_number
     */
    void index_wing_decomposition::decompose(const shared_ptr<abstract_bipartite_graph> &G,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                             const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                             double tau,
                                             const shared_ptr<thread_pool>& pool) {

        auto edge_set = G->get_edge_set();
        /**
         * @brief compute edge support
         */
        auto vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        vertex_priority_computation(G, vertex_priority_map, pool);
        edge_support_computation(G, edge_set,  edge_mutex_map, edge_support_map, vertex_priority_map, pool);

        auto B_mutex_map = make_shared<unordered_map<shared_ptr<priority_obeyed_bloom>,shared_ptr<mutex>>>();

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto rank_id = make_shared<uint32_t>(0);

        uint32_t k_max = estimate_maximal_wing(edge_support_map);
        uint32_t epsilon = k_max;
        while(epsilon > 0){
            auto subgraph = make_shared<abstract_bipartite_graph>();
            auto sub_edge_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
            auto sub_edge_support_map = make_shared<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>();

            for(const auto& e:*edge_set){
                if(edge_support_map->at(e) >= epsilon){
                    subgraph->insert_edge(e);
                    sub_edge_set->insert(e);
                    sub_edge_support_map->insert({e, 0});
                }
            }

            auto sub_vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();
            vertex_priority_computation(subgraph, sub_vertex_priority_map, pool);

            edge_support_computation(subgraph, sub_edge_set, edge_mutex_map, sub_edge_support_map, sub_vertex_priority_map, pool);
            remove_unsatisfied_edge(subgraph,sub_edge_set,sub_edge_support_map,epsilon);

            auto sub_bloom_index = compressed_index_construction(subgraph, sub_edge_set, sub_edge_support_map, edge_wing_map, edge_mutex_map, pool);
            sub_edge_support_map->clear();

            auto C = make_shared<unordered_map<shared_ptr<priority_obeyed_bloom>,uint32_t>>();
            {
                auto sub_bloom_map = sub_bloom_index->get_bloom_map();
                for(const auto&[p,B]:*sub_bloom_index->get_bloom_map())
                {
                    if(!B_mutex_map->count(B)){
                        B_mutex_map->insert({B, shared_ptr<mutex>()});
                    }
                }
                auto location_vector = pool->split_task(sub_bloom_map);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        auto &sub_begin = *location_vector->at(i);
                        auto &sub_end = *location_vector->at(i + 1);
                        for(auto iter = sub_begin; iter!=sub_end; ++iter){
                            auto &B = iter->second;
                            if(!B_mutex_map->at(B)){
                                B_mutex_map->at(B) = make_shared<mutex>();
                            }
                        }
                    });
                }
                for(const auto&[p,B]:*sub_bloom_index->get_bloom_map()){
                    C->insert({B,0});
                }
                pool->barrier();
            }

            while(!sub_edge_set->empty()){
                auto MBS = minimal_butterfly_support_finding(sub_bloom_index->get_edge_support_map(),edge_wing_map, pool);

                auto S_cur = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                scan(sub_edge_set, sub_bloom_index->get_edge_support_map(), edge_wing_map, MBS,S_cur, pool);

                while (!S_cur->empty()){
                    {
                        auto location_vector = pool->split_task(S_cur);
                        for(uint32_t i = 0; i< thread_number; ++i) {
                            pool->submit_task([=] {
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                    auto &e1 = *iter;
                                    edge_wing_map->at(e1) = sub_bloom_index->get_support(e1);
                                }
                            });
                        }
                        pool->barrier();
                    }

                    auto affected_B_set = make_shared<unordered_set<shared_ptr<priority_obeyed_bloom>>>();
                    {
                        set_edge_rank(S_cur, edge_rank_map, rank_id);
                        auto location_vector = pool->split_task(S_cur);
                        for(uint32_t i = 0; i < thread_number; ++i){
                            pool->submit_task([=]{
                                auto sub_affected_B_set = make_shared<unordered_set<shared_ptr<priority_obeyed_bloom>>>();

                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for(auto iter = sub_begin; iter != sub_end; ++iter){
                                    auto &e1 = *iter;

                                    edge_mutex_map->at(e1)->lock();
                                    auto B_set = sub_bloom_index->get_bloom_set(e1);

                                    for(const auto& B:*B_set){
                                        if(B->get_butterfly_count() > 0){
                                            B_mutex_map->at(B)->lock();
                                            auto e2 = B->get_twin(e1);
                                            B_mutex_map->at(B)->unlock();

                                            if(S_cur->count(e2) && edge_rank_map->at(e1) < edge_rank_map->at(e2)){
                                                continue;
                                            }

                                            sub_affected_B_set->insert(B);

                                            B_mutex_map->at(B)->lock();
                                            ++C->at(B);
                                            B_mutex_map->at(B)->unlock();

                                            auto k = B->get_k();
                                            if(edge_wing_map->at(e2) == 0){
                                                edge_mutex_map->at(e2)->lock();
                                                sub_bloom_index->update_support(e2,max(MBS,sub_bloom_index->get_support(e2)-(k-1)));
                                                sub_bloom_index->remove_edge(e2,B);
                                                edge_mutex_map->at(e2)->unlock();

                                                B_mutex_map->at(B)->lock();
                                                B->remove_edge(e2);
                                                B_mutex_map->at(B)->unlock();
                                            }
                                        }
                                    }
                                    edge_mutex_map->at(e1)->unlock();
                                }
                                global_mutex->lock();
                                affected_B_set->merge(*sub_affected_B_set);
                                global_mutex->unlock();
                            });
                        }
                        pool->barrier();
                    }

                    auto S_next = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                    {
                        auto location_vector = pool->split_task(affected_B_set);
                        for(uint32_t i = 0; i < thread_number; ++i) {
                            pool->submit_task([=] {
                                auto sub_S_next_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                                    auto &B = *iter;
                                    auto k = B->get_k();
                                    auto butterfly_count = (k - C->at(B)) * ((k - 1) - C->at(B)) / 2;
                                    B->set_butterfly_count(butterfly_count);
                                    if (butterfly_count == 0) {
                                        C->at(B) = k - 1;
                                    }

                                    for (const auto&[e1, e2]:*B->get_edge_map()) {
                                        if (!S_cur->count(e1)) {
                                            /**
                                             * @briref parallel visit
                                             */
                                            edge_mutex_map->at(e1)->lock();
                                            if (edge_wing_map->at(e1) == 0) {
                                                sub_bloom_index->update_support(e1, max(MBS,
                                                                                        sub_bloom_index->get_support(e1) -
                                                                                        C->at(B)));
                                            }
                                            edge_mutex_map->at(e1)->unlock();

                                            if (edge_wing_map->at(e1) == 0 && sub_bloom_index->get_support(e1) == MBS) {
                                                sub_S_next_set->insert(e1);
                                            }
                                        }
                                    }
                                }

                                global_mutex->lock();
                                S_next->merge(*sub_S_next_set);
                                global_mutex->unlock();
                            });
                        }
                        pool->barrier();
                    }

                    {
                        auto location_vector = pool->split_task(affected_B_set);
                        for(uint32_t i = 0; i < thread_number; ++i) {
                            pool->submit_task([=] {
                                auto &sub_begin = *location_vector->at(i);
                                auto &sub_end = *location_vector->at(i + 1);

                                for (auto iter = sub_begin; iter != sub_end; ++iter) {
                                    auto &B = *iter;
                                    C->at(B) = 0;
                                }
                            });
                        }
                        for(const auto&e:*S_cur){
                            sub_edge_set->erase(e);
                            sub_bloom_index->remove_edge(e);
                            edge_rank_map->at(e) = UINT32_MAX;
                        }
                        pool->barrier();
                    }
                    affected_B_set->clear();
                    S_cur->clear();
                    swap(*S_cur,*S_next);
                }
            }

            auto alpha = static_cast<uint32_t>(ceil((k_max) * tau));
            if(epsilon == 1)
            {
                epsilon = 0;
            }else
            {
                epsilon =  epsilon > alpha ? epsilon - alpha: 1;
            }
        }
    }

    /**
     * @details estimate the maximal wing number of a given graph
     * @param edge_wing_map
     * @return
     */
    uint32_t index_wing_decomposition::estimate_maximal_wing(
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map) {
        auto count_map = make_shared<map<uint32_t,uint32_t>>();
        for(const auto&[e,support]:*edge_wing_map){
            if(!count_map->count(support)){
                count_map->insert({support,0});
            }
            ++count_map->at(support);
        }
        uint32_t k_max = 0;
        uint32_t sum = 0;
        for(auto iter = count_map->rbegin();iter!=count_map->rend();++iter)
        {
            auto &[support,count] = *iter;
            sum += count;
            if(sum >= support){
                k_max = support;
                break;
            }
        }
        return k_max;
    }


    void index_wing_decomposition::init(const shared_ptr<abstract_bipartite_graph> &G,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        const shared_ptr<thread_pool>& pool)
    {
        auto edge_set = G->get_edge_set();
        for (const auto &e:*edge_set) {
            edge_mutex_map->insert({e, shared_ptr<mutex>()});
        }
        auto thread_number = pool->get_thread_number();
        auto location_vector = pool->split_task(edge_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &e = *iter;
                    edge_mutex_map->at(e) = make_shared<mutex>();
                }
            });
        }
        for (const auto &e:*edge_set) {
            edge_support_map->insert({e, 0});
            edge_rank_map->insert({e, UINT32_MAX});
            edge_wing_map->insert({e, 0});
        }
        pool->barrier();
    }

    uint32_t index_wing_decomposition::WS_init(const shared_ptr<abstract_bipartite_graph> &G,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                               const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                               const shared_ptr<thread_pool>& pool){
        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_rank_map->at(e) = UINT32_MAX;
            }
        });

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                WS->insert({e, 0});
            }
        });

        auto wing_edge_map = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
        for(const auto &[e, wing_number]:*edge_wing_map){
            if(!wing_edge_map->count(wing_number)){
                wing_edge_map->insert({wing_number, make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>()});
            }
            wing_edge_map->at(wing_number)->insert(e);
        }

        pool->barrier();

        auto thread_number = pool->get_thread_number();

        auto rank_id = make_shared<uint32_t>(0);
        for(const auto &[wing_number, e_set]:*wing_edge_map){
            set_edge_rank(e_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(e_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for(const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()){
                            if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
                                continue;
                            }

                            for(const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()){
                                if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if(!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                edge_mutex_map->at(e1)->lock();
                                ++WS->at(e1);
                                edge_mutex_map->at(e1)->unlock();

                                if(edge_wing_map->at(e2) == edge_wing_map->at(e1)){
                                    edge_mutex_map->at(e2)->lock();
                                    ++WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if(edge_wing_map->at(e3) == edge_wing_map->at(e1)){
                                    edge_mutex_map->at(e3)->lock();
                                    ++WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if(edge_wing_map->at(e4) == edge_wing_map->at(e1)){
                                    edge_mutex_map->at(e4)->lock();
                                    ++WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }
                            }
                        }
                    }
                });
            }
            pool->barrier();
        }
        return wing_edge_map->rbegin()->first;
    }

    uint32_t index_wing_decomposition::WS_WL_init(const shared_ptr<abstract_bipartite_graph> &G,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>,shared_ptr<mutex>>>& edge_mutex_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& edge_rank_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                  const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& WS,
                                                  const shared_ptr<unordered_map<uint32_t,shared_ptr<map<uint32_t,shared_ptr<unordered_set<uint32_t>>>>>> &WL,
                                                  const shared_ptr<thread_pool>& pool){

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                edge_rank_map->at(e) = UINT32_MAX;
            }
        });

        pool->submit_task([=]{
            for(const auto &[e, wing_number]:*edge_wing_map){
                WS->insert({e, 0});
            }
        });

        pool->submit_task([=]{
           WL->reserve(G->get_left_vertex_number()+G->get_right_vertex_number());
           for(const auto &[l, l_vertex]:*G->get_left_vertex_map()){
               WL->insert({l, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
           }
            for(const auto &[r, r_vertex]:*G->get_right_vertex_map()){
                WL->insert({r, shared_ptr<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>()});
            }
        });

        auto wing_edge_map = make_shared<map<uint32_t, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>>>();
        pool->submit_task([=]{
           for(const auto &[e, wing_number]:*edge_wing_map) {
               if(!wing_edge_map->count(wing_number)){
                   wing_edge_map->insert({wing_number, shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>()});
               }
           }
        });

        pool->barrier();

        auto thread_number = pool->get_thread_number();
        {
            {
                auto location_vector = pool->split_task(wing_edge_map);
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        for(auto iter = *location_vector->at(i); iter !=*location_vector->at(i + 1); ++iter){
                            pool->submit_task([=]{
                                auto wing_number = iter->first;
                                wing_edge_map->at(wing_number) = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();
                                for(const auto &[e, e_wing_number]:*edge_wing_map){
                                    if(e_wing_number == wing_number){
                                        wing_edge_map->at(wing_number)->insert(e);
                                    }
                                }
                            });
                        }
                    });
                }
                pool->barrier();
            }


            {
                auto location_vector = pool->split_task(G->get_left_vertex_map());
                for(uint32_t i = 0; i < thread_number; ++i){
                    pool->submit_task([=]{
                        for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                            auto &[l, l_vertex] = *iter;
                            WL->at(l) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                            for(const auto &[r, e]:*l_vertex->get_edge_map()){
                                auto wing_number = edge_wing_map->at(e);
                                if(!WL->at(l)->count(wing_number)){
                                    WL->at(l)->insert({wing_number, make_shared<unordered_set<uint32_t>>()});
                                }
                                WL->at(l)->at(wing_number)->insert(r);
                            }
                        }
                    });
                }
            }
            {
                auto location_vector = pool->split_task(G->get_right_vertex_map());
                for(uint32_t i = 0; i < pool->get_thread_number(); ++i){
                    pool->submit_task([=]{
                        for(auto iter = *location_vector->at(i); iter!=*location_vector->at(i + 1); ++iter){
                            auto &[r, r_vertex] = *iter;
                            WL->at(r) = make_shared<map<uint32_t, shared_ptr<unordered_set<uint32_t>>>>();
                            for(const auto &[l, e]:*r_vertex->get_edge_map()){
                                auto wing_number = edge_wing_map->at(e);
                                if(!WL->at(r)->count(wing_number)){
                                    WL->at(r)->insert({wing_number, make_shared<unordered_set<uint32_t>>()});
                                }
                                WL->at(r)->at(wing_number)->insert(l);
                            }
                        }
                    });
                }
            }
            pool->barrier();
        }



        auto rank_id = make_shared<uint32_t>(0);
        for(const auto &[wing_number, e_set]:*wing_edge_map){
            set_edge_rank(e_set, edge_rank_map, rank_id);
            auto location_vector = pool->split_task(e_set);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    for(auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter){
                        auto &e1 = *iter;

                        auto l1 = e1->get_left_vertex_id();
                        auto r1 = e1->get_right_vertex_id();

                        for(const auto &[r2, e2]:*G->get_left_vertex(l1)->get_edge_map()){
                            if(r2 == r1 || edge_rank_map->at(e2) < edge_rank_map->at(e1)){
                                continue;
                            }

                            for(const auto &[l2, e3]:*G->get_right_vertex(r1)->get_edge_map()){
                                if(l2 == l1 || edge_rank_map->at(e3) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                auto e4 = G->get_edge(l2, r2);
                                if(!e4 || edge_rank_map->at(e4) < edge_rank_map->at(e1)){
                                    continue;
                                }

                                edge_mutex_map->at(e1)->lock();
                                ++WS->at(e1);
                                edge_mutex_map->at(e1)->unlock();

                                if(edge_wing_map->at(e2) == edge_wing_map->at(e1)){
                                    edge_mutex_map->at(e2)->lock();
                                    ++WS->at(e2);
                                    edge_mutex_map->at(e2)->unlock();
                                }

                                if(edge_wing_map->at(e3) == edge_wing_map->at(e1)){
                                    edge_mutex_map->at(e3)->lock();
                                    ++WS->at(e3);
                                    edge_mutex_map->at(e3)->unlock();
                                }

                                if(edge_wing_map->at(e4) == edge_wing_map->at(e1)){
                                    edge_mutex_map->at(e4)->lock();
                                    ++WS->at(e4);
                                    edge_mutex_map->at(e4)->unlock();
                                }
                            }
                        }
                    }
                });
            }
            pool->barrier();
        }
        return wing_edge_map->rbegin()->first;
    }

    /**
     * @details construct bloom-index of left vertices
     * @param G
     * @param vertex_priority
     * @param edge_support_map
     * @param edge_mutex_map
     * @param bloom_index
     * @param bloom_index_mutex
     * @param thread_count
     */
    void index_wing_decomposition::left_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                           const shared_ptr<BE_index> &bloom_index,
                                                           const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto bloom_index_mutex = make_shared<mutex>();

        auto left_vertex_map = G->get_left_vertex_map();
        auto location_vector = pool->split_task(left_vertex_map);
        for (uint32_t i = 0; i < thread_number; ++i) {
            pool->submit_task([=] {
                auto sub_bloom_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                        shared_ptr<priority_obeyed_bloom>, hash_pair, equal_pair>>();

                for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                    auto &[l1, l1_vertex] = *iter;

                    auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    wedge_count_map->reserve(l1_vertex->get_edge_map()->size());

                    for (const auto &[r1, l1r1_edge]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority->at(r1) < vertex_priority->at(l1)) {
                            auto r1_vertex = G->get_right_vertex(r1);
                            for (const auto &[l2, l2r1_edge]: *r1_vertex->get_edge_map()) {
                                if (vertex_priority->at(l2) < vertex_priority->at(l1)) {
                                    if(!wedge_count_map->count(l2)){
                                        wedge_count_map->insert({l2, 0});
                                    }
                                    ++wedge_count_map->at(l2);
                                }
                            }
                        }
                    }

                    for (const auto &[r1, l1r1_edge]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority->at(r1) < vertex_priority->at(l1)) {
                            auto r1_vertex = G->get_right_vertex(r1);
                            for (const auto &[l2, l2r1_edge]: *r1_vertex->get_edge_map()) {
                                if (vertex_priority->at(l2) < vertex_priority->at(l1)) {
                                    if (wedge_count_map->at(l2) > 1) {
                                        if (!sub_bloom_map->count({l1, l2})) {
                                            auto B = make_shared<priority_obeyed_bloom>(l1, l2);

                                            auto butterfly_count =
                                                    wedge_count_map->at(l2) * (wedge_count_map->at(l2) - 1) / 2;
                                            B->set_butterfly_count(butterfly_count);

                                            sub_bloom_map->insert({{l1,l2}, B});
                                        }
                                        auto  B = sub_bloom_map->at({l1, l2});
                                        B->link_twin(l1r1_edge, l2r1_edge);

                                        edge_mutex_map->at(l1r1_edge)->lock();
                                        bloom_index->link_bloom(l1r1_edge, B);
                                        edge_mutex_map->at(l1r1_edge)->unlock();

                                        edge_mutex_map->at(l2r1_edge)->lock();
                                        bloom_index->link_bloom(l2r1_edge, B);
                                        edge_mutex_map->at(l2r1_edge)->unlock();
                                    }
                                }
                            }
                        }
                    }
                }

                bloom_index_mutex->lock();
                bloom_index->get_bloom_map()->merge(*sub_bloom_map);
                bloom_index_mutex->unlock();
            });
        }
        pool->barrier();
    }

    /**
     * @details construct sub-bloom-index of left vertices
     * @param G
     * @param vertex_priority
     * @param edge_support_map
     * @param edge_wing_map
     * @param edge_mutex_map
     * @param bloom_index
     * @param bloom_index_mutex
     * @param thread_count
     */
    void index_wing_decomposition::left_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                           const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                           const shared_ptr<BE_index> &bloom_index,
                                                           const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto left_vertex_map = G->get_left_vertex_map();
        auto location_vector = pool->split_task(left_vertex_map);

        for (uint32_t i = 0; i < thread_number; ++i) {
            pool->submit_task([=] {
                auto sub_bloom_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                        shared_ptr<priority_obeyed_bloom>, hash_pair, equal_pair>>();

                for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                    auto &[l1, l1_vertex] = *iter;
                    auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    for (const auto &[r1, e1]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority->at(r1) < vertex_priority->at(l1)) {
                            auto r1_vertex = G->get_right_vertex(r1);
                            for (const auto &[l2, e2]: *r1_vertex->get_edge_map()) {
                                if (vertex_priority->at(l2) < vertex_priority->at(l1)) {
                                    if (!wedge_count_map->count(l2)) {
                                        wedge_count_map->insert({l2, 0});
                                    }
                                    ++wedge_count_map->at(l2);
                                }
                            }
                        }
                    }
                    for (const auto &[r1, e1]: *l1_vertex->get_edge_map()) {
                        if (vertex_priority->at(r1) < vertex_priority->at(l1)) {
                            auto r1_vertex = G->get_right_vertex(r1);
                            for (const auto &[l2, e2]: *r1_vertex->get_edge_map()) {
                                if (vertex_priority->at(l2) < vertex_priority->at(l1)) {
                                    if (wedge_count_map->at(l2) > 1) {

                                        if (!sub_bloom_map->count({l1, l2})) {
                                            auto B = make_shared<priority_obeyed_bloom>(l1, l2);
                                            auto butterfly_count =
                                                    wedge_count_map->at(l2) * (wedge_count_map->at(l2) - 1) / 2;
                                            B->set_butterfly_count(butterfly_count);

                                            sub_bloom_map->insert({{l1, l2}, B});
                                        }
                                        auto B = sub_bloom_map->at({l1, l2});
                                        B->link_twin(e1, e2);

                                        if (edge_wing_map->at(e1) == 0) {
                                            edge_mutex_map->at(e1)->lock();
                                            bloom_index->link_bloom(e1, B);
                                            edge_mutex_map->at(e1)->unlock();
                                        }

                                        if (edge_wing_map->at(e2) == 0) {

                                            edge_mutex_map->at(e2)->lock();
                                            bloom_index->link_bloom(e2, B);
                                            edge_mutex_map->at(e2)->unlock();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                global_mutex->lock();
                bloom_index->get_bloom_map()->merge(*sub_bloom_map);
                global_mutex->unlock();
            });
        }
        pool->barrier();

    }

    /**
     * @details find the minimal support of edges
     * @param edge_support_map
     * @param edge_wing_map
     * @return
     */
    uint32_t index_wing_decomposition::minimal_butterfly_support_finding(
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map) {
        auto MBS = edge_support_map->begin()->second;
        for(const auto&[e,support]:*edge_support_map){
            if(edge_wing_map->at(e) > 0){
                continue;
            }
            if(support < MBS)
            {
                MBS = support;
            }
        }
        return MBS;
    }

    uint32_t index_wing_decomposition::minimal_butterfly_support_finding(
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
            const shared_ptr<thread_pool>& pool) {
        auto MBS = make_shared<uint32_t>(edge_support_map->begin()->second);
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto location_vector = pool->split_task(edge_support_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                uint32_t sub_MBS =UINT32_MAX;

                for(auto iter = sub_begin; iter!=sub_end; ++iter){
                    auto &e = iter->first;
                    if(edge_wing_map->at(e) > 0){
                        continue;
                    }
                    if(edge_support_map->at(e) < sub_MBS)
                    {
                        sub_MBS = edge_support_map->at(e);
                    }
                }

                global_mutex->lock();
                if(sub_MBS < *MBS){
                    *MBS = sub_MBS;
                }
                global_mutex->unlock();
            });
        }
        pool->barrier();

        return *MBS;
    }


    /**
     * @details remove unsatisfied edges in a epsilon-wing
     * @param G
     * @param edge_set
     * @param edge_support_map
     * @param epsilon
     */
    void index_wing_decomposition::remove_unsatisfied_edge(const shared_ptr<abstract_bipartite_graph> &G,
                                                           const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& edge_set,
                                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                           uint32_t epsilon)
    {
        auto evict_set = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

        for(const auto&e:*edge_set){
            if(edge_support_map->at(e) < epsilon)
            {
                evict_set->insert(e);
            }
        }
        while (!evict_set->empty())
        {
            auto e1 = *evict_set->begin();
            evict_set->erase(e1);

            auto l1 = e1->get_left_vertex_id();
            auto r1 = e1->get_right_vertex_id();
            auto l1_vertex = G->get_left_vertex(l1);
            auto r1_vertex = G->get_right_vertex(r1);
            for(const auto&[r2,e2]:*l1_vertex->get_edge_map())
            {
                if(r2 == r1){
                    continue;
                }
                for(const auto&[l2,e3]:*r1_vertex->get_edge_map()){
                    if(l2 == l1){
                        continue;
                    }
                    auto e4 = G->get_edge(l2,r2);
                    if(e4){
                        if(edge_support_map->at(e2) >= epsilon){
                            --edge_support_map->at(e2);
                            if(edge_support_map->at(e2) < epsilon)
                            {
                                evict_set->insert(e2);
                            }
                        }

                        if(edge_support_map->at(e3) >= epsilon)
                        {
                            --edge_support_map->at(e3);
                            if(edge_support_map->at(e3) < epsilon)
                            {
                                evict_set->insert(e3);
                            }
                        }

                        if(edge_support_map->at(e4) >= epsilon)
                        {
                            --edge_support_map->at(e4);
                            if(edge_support_map->at(e4) < epsilon)
                            {
                                evict_set->insert(e4);
                            }
                        }
                    }
                }
            }
            G->remove_edge(e1);
            edge_set->erase(e1);
            edge_support_map->erase(e1);
        }
    }


    /**
     * @details construct the bloom-index of right vertices
     * @param G
     * @param vertex_priority_map
     * @param edge_support_map
     * @param edge_mutex_map
     * @param bloom_index
     * @param global_mutex
     * @param thread_count
     */
    void index_wing_decomposition::right_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                            const shared_ptr<BE_index> &bloom_index,
                                                            const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto right_vertex_map = G->get_right_vertex_map();
        auto location_vector = pool->split_task(right_vertex_map);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_bloom_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                        shared_ptr<priority_obeyed_bloom>, hash_pair, equal_pair>>();

                for(auto iter = *location_vector->at(i); iter!= *location_vector->at(i + 1); ++iter){
                    auto &[r1,r1_vertex] = *iter;

                    auto wedge_count_map = make_shared<unordered_map<uint32_t, uint32_t>>();
                    wedge_count_map->reserve(r1_vertex->get_edge_map()->size());

                    for(const auto&[l1,l1r1_edge]:*r1_vertex->get_edge_map()){
                        if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)){
                            auto l1_vertex = G->get_left_vertex(l1);
                            for(const auto&[r2,l1r2_edge]:*l1_vertex->get_edge_map()){
                                if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1))
                                {
                                    if(!wedge_count_map->count(r2))
                                    {
                                        wedge_count_map->insert({r2,0});
                                    }
                                    ++wedge_count_map->at(r2);
                                }
                            }
                        }
                    }
                    for(const auto&[l1,e1]:*r1_vertex->get_edge_map()){
                        if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)){
                            auto l1_vertex = G->get_left_vertex(l1);
                            for(const auto&[r2,e2]:*l1_vertex->get_edge_map()){
                                if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1)){
                                    if(wedge_count_map->at(r2) > 1){
                                        if(!sub_bloom_map->count({r1, r2})){
                                            auto B = make_shared<priority_obeyed_bloom>(r1,r2);
                                            sub_bloom_map->insert({{r1, r2},B});
                                        }
                                        auto B = sub_bloom_map->at({r1,r2});

                                        auto butterfly_count = wedge_count_map->at(r2)*(wedge_count_map->at(r2)-1)/2;
                                        B->set_butterfly_count(butterfly_count);
                                        B->link_twin(e1, e2);

                                        edge_mutex_map->at(e1)->lock();
                                        bloom_index->link_bloom(e1, B);
                                        edge_mutex_map->at(e1)->unlock();

                                        edge_mutex_map->at(e2)->lock();
                                        bloom_index->link_bloom(e2, B);
                                        edge_mutex_map->at(e2)->unlock();
                                    }
                                }
                            }
                        }
                    }
                }

                global_mutex->lock();
                bloom_index->get_bloom_map()->merge(*sub_bloom_map);
                global_mutex->unlock();
            });
        }
        pool->barrier();
    }

    /**
     * @details construct the sub-bloom-index of right vertices
     * @param G
     * @param vertex_priority
     * @param edge_support_map
     * @param edge_wing_map
     * @param edge_mutex_map
     * @param bloom_index
     * @param bloom_index_mutex
     * @param thread_count
     */
    void index_wing_decomposition::right_index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                            const shared_ptr<BE_index> &bloom_index,
                                                            const shared_ptr<thread_pool>& pool) {
        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto right_vertex_map = G->get_right_vertex_map();
        auto location_vector = pool->split_task(right_vertex_map);
        for(uint32_t i = 0; i < thread_number; ++i) {
            pool->submit_task([=] {
                auto sub_bloom_map = make_shared<unordered_map<pair<uint32_t, uint32_t>,
                        shared_ptr<priority_obeyed_bloom>, hash_pair, equal_pair>>();

                for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                    auto &[r1, r1_vertex] = *iter;

                    auto wedge_count_map = make_shared<unordered_map<uint32_t,uint32_t>>();
                    for(const auto&[l1,l1r1_edge]:*r1_vertex->get_edge_map()){
                        if(vertex_priority->at(l1) < vertex_priority->at(r1)){
                            auto l1_vertex = G->get_left_vertex(l1);
                            for(const auto&[r2,l1r2_edge]:*l1_vertex->get_edge_map()){
                                if(vertex_priority->at(r2) < vertex_priority->at(r1))
                                {
                                    if(!wedge_count_map->count(r2))
                                    {
                                        wedge_count_map->insert({r2,0});
                                    }
                                    ++wedge_count_map->at(r2);
                                }
                            }
                        }
                    }
                    for(const auto&[l1,l1r1_edge]:*r1_vertex->get_edge_map()){
                        if(vertex_priority->at(l1) < vertex_priority->at(r1)){
                            auto l1_vertex = G->get_left_vertex(l1);
                            for(const auto&[r2,l1r2_edge]:*l1_vertex->get_edge_map()){
                                if(vertex_priority->at(r2) < vertex_priority->at(r1)){
                                    if(wedge_count_map->at(r2) > 1){

                                        if(!sub_bloom_map->count({r1,r2})){
                                            auto B = make_shared<priority_obeyed_bloom>(r1,r2);
                                            sub_bloom_map->insert({{r1,r2}, B});
                                        }
                                        auto B = sub_bloom_map->at({r1,r2});

                                        auto butterfly_count = wedge_count_map->at(r2)*(wedge_count_map->at(r2)-1)/2;
                                        B->set_butterfly_count(butterfly_count);
                                        B->link_twin(l1r1_edge, l1r2_edge);

                                        if(edge_wing_map->at(l1r1_edge) == 0)
                                        {
                                            edge_mutex_map->at(l1r1_edge)->lock();
                                            bloom_index->link_bloom(l1r1_edge, B);
                                            edge_mutex_map->at(l1r1_edge)->unlock();
                                        }

                                        if(edge_wing_map->at(l1r2_edge) == 0)
                                        {
                                            edge_mutex_map->at(l1r2_edge)->lock();
                                            bloom_index->link_bloom(l1r2_edge, B);
                                            edge_mutex_map->at(l1r2_edge)->unlock();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                global_mutex->lock();
                bloom_index->get_bloom_map()->merge(*sub_bloom_map);
                global_mutex->unlock();
            });
        }
        pool->barrier();
    }

    /**
     * @details construct the bloom-index of a given graph
     * @param G
     * @param edge_set
     * @param edge_wing_map
     * @param edge_rank_map
     * @param edge_mutex_map
     * @param bloom_index_mutex
     * @param thread_number
     * @return
     */
    shared_ptr<BE_index> index_wing_decomposition::index_construction(const shared_ptr<abstract_bipartite_graph> &G,
                                                                      const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                                                      const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                                      const shared_ptr<thread_pool>& pool) {

        auto bloom_index = make_shared<BE_index>();
        auto vertex_priority_map = make_shared<unordered_map<uint32_t, uint32_t>>();
        vertex_priority_computation(G, vertex_priority_map, pool);

        edge_support_computation(G, edge_set, edge_mutex_map,  edge_support_map, vertex_priority_map, pool);

        for(auto iter = edge_set->begin();iter!= edge_set->end();){
            auto e = *iter;
            ++iter;
            edge_wing_map->insert({e,0});
            auto support = edge_support_map->at(e);
            if(support > 0){
                bloom_index->insert_edge(e,support);
            }else
            {
                edge_set->erase(e);
                edge_support_map->erase(e);
            }
        }

        left_index_construction(G, vertex_priority_map, edge_support_map, edge_mutex_map, bloom_index, pool);
        right_index_construction(G, vertex_priority_map, edge_support_map, edge_mutex_map, bloom_index, pool);
        return bloom_index;
    }


    /**
     * @details find a a set of edges whose support equal to the minimal support
     * @param edge_set
     * @param edge_support_map
     * @param edge_wing_map
     * @param MBS
     * @param thread_number
     */
    void index_wing_decomposition::scan(const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        uint32_t MBS,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& S_cur,
                                        const shared_ptr<thread_pool>& pool)
    {

        auto thread_number = pool->get_thread_number();
        auto global_mutex = make_shared<mutex>();

        auto location_vector = pool->split_task(edge_set);
        for(uint32_t i = 0; i < thread_number; ++i){
            pool->submit_task([=]{
                auto sub_S_cur = make_shared<unordered_set<shared_ptr<abstract_bipartite_edge>>>();

                auto &sub_begin = *location_vector->at(i);
                auto &sub_end = *location_vector->at(i + 1);

                for(auto iter = sub_begin; iter !=sub_end; ++iter){
                    auto &e = *iter;
                    if (edge_wing_map->at(e) == 0 && edge_support_map->at(e) == MBS) {
                        sub_S_cur->insert(e);
                    }
                }

                global_mutex->lock();
                S_cur->merge(*sub_S_cur);
                global_mutex->unlock();
            });
        }
        pool->barrier();
    }

    void index_wing_decomposition::scan(const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_wing_map,
                                        uint32_t MBS,
                                        const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>>& S_cur)
    {
        for(const auto&e:*edge_set){
            if (edge_wing_map->at(e) == 0 && edge_support_map->at(e) == MBS) {
                S_cur->insert(e);
            }
        }
    }

    void index_wing_decomposition::set_edge_rank(const shared_ptr<unordered_set<shared_ptr<scnu::abstract_bipartite_edge>>> &edge_set,
                                                 const shared_ptr<unordered_map<shared_ptr<scnu::abstract_bipartite_edge>, uint32_t>> &edge_rank_map,
                                                 const shared_ptr<uint32_t>& rank_id) {
        for(const auto &e:*edge_set){
            *rank_id = *rank_id + 1;
            edge_rank_map->at(e) = *rank_id;
        }
    }


    void index_wing_decomposition::edge_support_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                            const shared_ptr<unordered_set<shared_ptr<abstract_bipartite_edge>>> &edge_set,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, shared_ptr<mutex>>> &edge_mutex_map,
                                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>> &edge_support_map,
                                                            const shared_ptr<unordered_map<uint32_t, uint32_t>> &vertex_priority_map,
                                                            const shared_ptr<thread_pool> &pool) {
        auto thread_number = pool->get_thread_number();
        {
            /**
             * @brief compute wedges from left vertices
             */
            auto left_vertex_map = G->get_left_vertex_map();
            auto location_vector = pool->split_task(left_vertex_map);

            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter !=sub_end; ++iter){
                        auto& [l1,l1_vertex] = *iter;
                        auto wedge_count_map = make_shared<unordered_map<uint32_t,uint32_t>>();

                        for(const auto&[r1,e1]:*l1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(r1) < vertex_priority_map->at(l1)){
                                auto r1_vertex = G->get_right_vertex(r1);
                                for(const auto&[l2,e2]:*r1_vertex->get_edge_map()){
                                    if(vertex_priority_map->at(l2) < vertex_priority_map->at(l1)){
                                        if(!wedge_count_map->count(l2)){
                                            wedge_count_map->insert({l2, 0});
                                        }
                                        ++wedge_count_map->at(l2);
                                    }
                                }
                            }
                        }
                        for(const auto&[r1,e1_edge]:*l1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(r1) < vertex_priority_map->at(l1)) {
                                auto r1_vertex = G->get_right_vertex(r1);
                                for (const auto&[l2, e2]:*r1_vertex->get_edge_map()) {
                                    if(vertex_priority_map->at(l2) < vertex_priority_map->at(l1)){
                                        if(wedge_count_map->at(l2) > 1){
                                            auto delta = wedge_count_map->at(l2) - 1;

                                            edge_mutex_map->at(e1_edge)->lock();
                                            edge_support_map->at(e1_edge) += delta;
                                            edge_mutex_map->at(e1_edge)->unlock();

                                            edge_mutex_map->at(e2)->lock();
                                            edge_support_map->at(e2) += delta;
                                            edge_mutex_map->at(e2)->unlock();
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            }
        }

        {
            /**
             * @brief compute wedges from right vertices
             */
            auto right_vertex_map = G->get_right_vertex_map();
            auto location_vector = pool->split_task(right_vertex_map);

            for(uint32_t i = 0; i < thread_number; ++i) {
                pool->submit_task([=] {
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = *location_vector->at(i + 1);

                    for (auto iter = sub_begin; iter != sub_end; ++iter) {
                        auto& [r1, r1_vertex] = *iter;
                        auto wedge_count_map = make_shared<unordered_map<uint32_t,uint32_t>>();

                        for(const auto&[l1, e1]:*r1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)){
                                auto l1_vertex = G->get_left_vertex(l1);
                                for(const auto&[r2,e2]:*l1_vertex->get_edge_map()){
                                    if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1)){
                                        if(!wedge_count_map->count(r2)){
                                            wedge_count_map->insert({r2, 0});
                                        }
                                        ++wedge_count_map->at(r2);
                                    }
                                }
                            }
                        }
                        for(const auto&[l1,e1]:*r1_vertex->get_edge_map()){
                            if(vertex_priority_map->at(l1) < vertex_priority_map->at(r1)) {
                                auto l1_vertex = G->get_left_vertex(l1);
                                for (const auto&[r2, e2]:*l1_vertex->get_edge_map()) {
                                    if(vertex_priority_map->at(r2) < vertex_priority_map->at(r1)){
                                        if(wedge_count_map->at(r2) > 1){
                                            auto delta = wedge_count_map->at(r2)-1;

                                            edge_mutex_map->at(e1)->lock();
                                            edge_support_map->at(e1) += delta;
                                            edge_mutex_map->at(e1)->unlock();

                                            edge_mutex_map->at(e2)->lock();
                                            edge_support_map->at(e2) += delta;
                                            edge_mutex_map->at(e2)->unlock();
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            }
        }
        pool->barrier();
    }


    /**
     * @details assign_value a priority for each vertex
     * @note each vertex must have a unique identifier
     * @param G
     * @return
     */
    void index_wing_decomposition::vertex_priority_computation(const shared_ptr<abstract_bipartite_graph> &G,
                                                               const shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_priority_map,
                                                               const shared_ptr<thread_pool> &pool) {
        auto global_mutex = make_shared<mutex>();
        auto max_degree = make_shared<uint32_t>(0);
        {
            {
                auto left_vertex_map = G->get_left_vertex_map();
                auto location_vector = pool->split_task(left_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        uint32_t sub_max = 0;
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [l, l_vertex] = *iter;
                            if (l_vertex->get_degree() > sub_max) {
                                sub_max = l_vertex->get_degree();
                            }
                        }

                        global_mutex->lock();
                        if (sub_max > *max_degree) {
                            *max_degree = sub_max;
                        }
                        global_mutex->unlock();
                    });
                }
            }

            {
                auto right_vertex_map = G->get_right_vertex_map();
                auto location_vector = pool->split_task(right_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        uint32_t sub_max = 0;
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [r, r_vertex] = *iter;
                            if (r_vertex->get_degree() > sub_max) {
                                sub_max = r_vertex->get_degree();
                            }
                        }

                        global_mutex->lock();
                        if (sub_max > *max_degree) {
                            *max_degree = sub_max;
                        }
                        global_mutex->unlock();
                    });
                }
            }
            pool->barrier();
        }

        auto degree_vector = make_shared<vector<shared_ptr<set<uint32_t>>>>(*max_degree + 1);
        auto degree_mutex_vector = make_shared<vector<shared_ptr<mutex>>>(*max_degree + 1);
        {
            auto location_vector = pool->split_task(degree_mutex_vector);
            for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                pool->submit_task([=] {
                    for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                        *iter = make_shared<mutex>();
                    }
                });
            }
            pool->barrier();
        }
        {
            {
                auto left_vertex_map = G->get_left_vertex_map();
                auto location_vector = pool->split_task(left_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [l, l_vertex] = *iter;
                            auto degree = l_vertex->get_degree();
                            degree_mutex_vector->at(degree)->lock();
                            if (!degree_vector->at(degree)) {
                                degree_vector->at(degree) = make_shared<set<uint32_t>>();
                            }
                            degree_vector->at(degree)->insert(l);
                            degree_mutex_vector->at(degree)->unlock();
                        }
                    });
                }
            }

            {
                auto right_vertex_map = G->get_right_vertex_map();
                auto location_vector = pool->split_task(right_vertex_map);
                for (uint32_t i = 0; i < pool->get_thread_number(); ++i) {
                    pool->submit_task([=] {
                        for (auto iter = *location_vector->at(i); iter != *location_vector->at(i + 1); ++iter) {
                            auto [r, r_vertex] = *iter;
                            auto degree = r_vertex->get_degree();
                            degree_mutex_vector->at(degree)->lock();
                            if (!degree_vector->at(degree)) {
                                degree_vector->at(degree) = make_shared<set<uint32_t>>();
                            }
                            degree_vector->at(degree)->insert(r);
                            degree_mutex_vector->at(degree)->unlock();
                        }
                    });
                }
            }
            pool->barrier();
        }
        uint32_t i = 0;
        vertex_priority_map->reserve(G->get_left_vertex_number() + G->get_right_vertex_number());
        for (auto iter1 = degree_vector->rbegin(); iter1 != degree_vector->rend(); ++iter1) {
            if (*iter1) {
                auto v_set = *iter1;
                for (auto iter2 = v_set->rbegin(); iter2 != v_set->rend(); ++iter2) {
                    vertex_priority_map->insert({*iter2, ++i});
                }
            }
        }
    }
}
