
#include "multiple_core/multiple_core.h"

namespace scnu{
    multiple_core::multiple_core(
            const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> &other_vertex_index_map) :
            basic_vertex_index_map(other_vertex_index_map){

    }

    shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>> multiple_core::get_core_pair_set(uint32_t delta) {
        auto core_pair_map = std::make_shared<unordered_set<pair<uint32_t, uint32_t>, hash_pair>>();
        for(const auto &[u, u_index]:*basic_vertex_index_map){
            for(const auto &[h, k]:*u_index->get_left_map()){
                if(h < delta || k < delta){
                    continue;
                }
                for(uint32_t index = h; index <= k; ++index){
                    core_pair_map->insert({index, h});
                }
            }

            for(const auto &[k, h]:*u_index->get_right_map()){
                if(k < delta || h < delta){
                    continue;
                }
                for(uint32_t index = k; index <= h; ++index){
                    core_pair_map->insert({k, index});
                }
            }
        }
        return core_pair_map;
    }

    uint32_t multiple_core::get_maximal_delta() {
        uint32_t max_delta = 0;
        for(const auto &[u, u_index]:*basic_vertex_index_map){
            auto delta = u_index->get_delta();
            if(delta > max_delta){
                max_delta = delta;
            }
        }
        return max_delta;
    }

    uint32_t multiple_core::get_maximal_h() {
        uint32_t max_h = 0;
        for(const auto &[u, u_index]:*basic_vertex_index_map){
            auto h = u_index->get_h(1);
            if(h > max_h){
                max_h = h;
            }
        }
        return max_h;
    }

    uint32_t multiple_core::get_maximal_h(uint32_t k) {
        uint32_t max_h = 0;
        for(const auto &[u, u_index]:*basic_vertex_index_map){
            auto h = u_index->get_h(k);
            if(h > max_h){
                max_h = h;
            }
        }
        return max_h;
    }

    uint32_t multiple_core::get_maximal_k() {
        uint32_t max_k = 0;
        for(const auto &[u, u_index]:*basic_vertex_index_map){
            auto k = u_index->get_k(1);
            if(k > max_k){
                max_k = k;
            }
        }
        return max_k;
    }

    uint32_t multiple_core::get_maximal_k(uint32_t h) {
        uint32_t max_k = 0;
        for(const auto &[u, u_index]:*basic_vertex_index_map){
            auto k = u_index->get_k(h);
            if(k > max_k){
                max_k = k;
            }
        }
        return max_k;
    }
}