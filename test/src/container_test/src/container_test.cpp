
#include "container_test/container_test.h"


int main()
{

    struct pair_compare{
        bool operator()(const pair<uint32_t,uint32_t>& p1,const pair<uint32_t,uint32_t>& p2) const{
            return p1.second < p2.second;
        }
    };

    auto vertex_degree_map = make_shared<unordered_map<uint32_t, uint32_t>>();
    for(uint32_t i = 1; i <= 10;++i){
        vertex_degree_map->insert({i, 10-i});
    }

    auto degree_pair_set = make_shared<std::multiset<pair<uint32_t,uint32_t>, pair_compare>>();
    for(const auto &[v,v_degree]:*vertex_degree_map){
        degree_pair_set->insert({v,v_degree});
    }

    auto p = std::make_pair(5, vertex_degree_map->at(5));
    degree_pair_set->erase(p);
    p.second += 3;
    degree_pair_set->insert(p);

    for(const auto&p:*degree_pair_set){
        printf("%u,%u\n",p.first,p.second);
    }

}