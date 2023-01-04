
#include "bipartite_core/bipartite_core_compare.h"

namespace scnu{
    void bipartite_core_compare::convert(
            const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
            const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map)
    {
        for(const auto&[p,core]:*core_map){
            auto [i,j] = p;
            for(const auto&l:*core->get_left_vertex_set()){
                if(!left_index_map->count(l)){
                    left_index_map->insert({l,make_shared<left_vertex_index>()});
                }
                left_index_map->at(l)->insert(i, j);
            }

            for(const auto&r:*core->get_right_vertex_set()){
                if(!right_index_map->count(r)){
                    right_index_map->insert({r,make_shared<right_vertex_index>()});
                }
                right_index_map->at(r)->insert(j, i);
            }
        }
    }

    void bipartite_core_compare::convert(const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                        const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                        const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map)
    {
        for(const auto &[l,l_node]:*left_index_map)
        {
            auto index_map = l_node->get_index_map();
            for(uint32_t i = 1;i<= index_map->size();++i){
                for(uint32_t j = 1;j<=index_map->at(i);++j){
                    if(!core_map->count({i,j})){
                        core_map->insert({make_pair(i,j),make_shared<bipartite_core>()});
                    }
                    core_map->at({i,j})->insert_left_vertex(l);
                }
            }
        }

        for(const auto& [r,r_node]:*right_index_map){
            auto index_map = r_node->get_index_map();
            for(uint32_t j = 1; j<= index_map->size();++j){
                for(uint32_t i =1;i<=index_map->at(j);++i){
                    if(!core_map->count({i,j})){
                        core_map->insert({make_pair(i,j),make_shared<bipartite_core>()});
                    }
                    core_map->at({i,j})->insert_right_vertex(r);
                }
            }
        }
    }

    void bipartite_core_compare::convert(const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                                         const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map,
                                         const shared_ptr<unordered_set<pair<uint32_t, uint32_t>,hash_pair, equal_pair>> &core_set)
    {
        for(const auto &[l,l_node]:*left_index_map)
        {
            auto index_map = l_node->get_index_map();
            for(uint32_t i = 1;i<= index_map->size();++i){
                for(uint32_t j = 1;j<=index_map->at(i);++j){
                    if(!core_set->count({i,j})){
                        core_set->insert({make_pair(i,j)});
                    }
                }
            }
        }

        for(const auto& [r,r_node]:*right_index_map){
            auto index_map = r_node->get_index_map();
            for(uint32_t j = 1; j<= index_map->size();++j){
                for(uint32_t i =1;i<=index_map->at(j);++i){
                    if(!core_set->count({i,j})){
                        core_set->insert({i,j});
                    }
                }
            }
        }
    }


    bool bipartite_core_compare::sub(const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map1,
                    const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map2) {
        for(const auto&[p,core]:*core_map1){
            if(!core_map2->count(p)){
                printf("Unequal size!\n");
                return false;
            }
            for(const auto&l:*core->get_left_vertex_set()){
                if(!core_map2->at(p)->count_left_vertex(l)){
                    printf("%u,%u,%u\n",l,p.first,p.second);
                    return false;
                }
            }

            for(const  auto&r:*core->get_right_vertex_set()){
                if(!core_map2->at(p)->count_right_vertex(r)){
                    printf("%u,%u,%u",r,p.first,p.second);
                    return false;
                }
            }
        }
        return true;
    }

    bool bipartite_core_compare::same(const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map1,
                     const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map2){
        return sub(core_map1,core_map2) && sub(core_map2,core_map1);
    }

    bool bipartite_core_compare::same(const shared_ptr<unordered_map<pair<uint32_t, uint32_t>, shared_ptr<bipartite_core>, hash_pair, equal_pair>> &core_map,
                     const shared_ptr<unordered_map<uint32_t, shared_ptr<left_vertex_index>>> &left_index_map,
                     const shared_ptr<unordered_map<uint32_t, shared_ptr<right_vertex_index>>> &right_index_map)
    {
        auto contrastive_left_index_map = make_shared<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>();
        auto contrastive_right_index_map = make_shared<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>();
        convert(core_map, contrastive_left_index_map, contrastive_right_index_map);
        return same(contrastive_left_index_map, contrastive_right_index_map, left_index_map, right_index_map);
    }


    bool bipartite_core_compare::sub(const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map1,
             const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map1,
             const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map2,
             const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map2) {

        for(const auto& [l,l_node1]:*left_index_map1)
        {
            if(!left_index_map2->count(l)){
                printf("\nUnequal Left Key: %u\n",l);
                return false;
            }
            auto l_node2 = left_index_map2->at(l);
            if(!l_node1->compare(l_node2))
            {
                printf("\nUnequal Left Index: %u!\n", l);
                return false;
            }
        }

        for(const auto& [r,r_node1]:*right_index_map1)
        {
            if(!right_index_map2->count(r))
            {
                printf("\nUnequal Right Key: %u\n",r);
                return false;
            }
            auto r_node2 = right_index_map2->at(r);
            if(!r_node1->compare(r_node2)){
                printf("\nUnequal Right Index: %u!\n", r);
                return false;
            }
        }
        return true;
    }

    bool bipartite_core_compare::same(const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map1,
              const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map1,
              const shared_ptr<unordered_map<uint32_t,shared_ptr<left_vertex_index>>>& left_index_map2,
              const shared_ptr<unordered_map<uint32_t,shared_ptr<right_vertex_index>>>& right_index_map2){
        return sub(left_index_map1,right_index_map1,left_index_map2,right_index_map2)
               && sub(left_index_map2,right_index_map2, left_index_map1, right_index_map1);
    }
}
