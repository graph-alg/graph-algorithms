/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : temporal_core_utility.h
* @brief      : a graph structure for temporal core
* @version    : 1.0
* @date       : 2021/08/06
******************************************************************************************************************/

#pragma once
#include "multiple_core/multiple_core_pair_map_index.h"
#include "multiple_core/multiple_core_pair_set_index.h"
#include "multiple_core/multiple_core_number_set_index.h"

namespace scnu{
    class multiple_core {
    public:
        multiple_core(const shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>>& other_vertex_index_map);

        template<class index_container_type>
        void convert(const shared_ptr<unordered_map<uint32_t, shared_ptr<index_container_type>>>& destination_vertex_index_map){
            destination_vertex_index_map->reserve(basic_vertex_index_map->size());
            for(const auto &[u,u_index]:*basic_vertex_index_map){
                destination_vertex_index_map->insert({u, make_shared<index_container_type>()});
            }

            for(const auto &[u,u_index]:*basic_vertex_index_map){
                for(const auto &[h, k]:*u_index->get_left_map()){
                    for(uint32_t index = h; index <= k; ++index){
                        destination_vertex_index_map->at(u)->insert(index, h);
                    }
                }

                for(const auto &[k, h]:*u_index->get_right_map()){
                    for(uint32_t index = k; index <= h; ++index){
                        destination_vertex_index_map->at(u)->insert(k, index);
                    }
                }
            }
        }

        template<class index_container_type>
        void convert(const shared_ptr<unordered_map<uint32_t, shared_ptr<index_container_type>>>& destination_vertex_index_map,
                     const shared_ptr<thread_pool>& pool){
            destination_vertex_index_map->reserve(basic_vertex_index_map->size());
            for(const auto &[u,u_index]:*basic_vertex_index_map){
                destination_vertex_index_map->insert({u, make_shared<index_container_type>()});
            }

            auto thread_number = pool->get_thread_number();
            auto location_vector = pool->split_task(basic_vertex_index_map);
            for(uint32_t i = 0; i < thread_number; ++i){
                pool->submit_task([=]{
                    auto &sub_begin = *location_vector->at(i);
                    auto &sub_end = * location_vector->at(i + 1);

                    for(auto iter = sub_begin; iter!=sub_end; ++iter){
                        auto &[u, u_index] = *iter;

                        for(const auto &[h, k]:*u_index->get_left_map()){
                            for(uint32_t index = h; index <= k; ++index){
                                destination_vertex_index_map->at(u)->insert(index, h);
                            }
                        }

                        for(const auto &[k, h]:*u_index->get_right_map()){
                            for(uint32_t index = k; index <= h; ++index){
                                destination_vertex_index_map->at(u)->insert(k, index);
                            }
                        }
                    }
                });
            }
            pool->barrier();
        }


        shared_ptr<unordered_set<pair<uint32_t, uint32_t>, hash_pair>> get_core_pair_set(uint32_t delta);

        uint32_t get_maximal_delta();

        uint32_t  get_maximal_h();

        uint32_t  get_maximal_h(uint32_t k);

        uint32_t  get_maximal_k();

        uint32_t  get_maximal_k(uint32_t h);



    private:
        shared_ptr<unordered_map<uint32_t, shared_ptr<multiple_core_pair_map_index>>> basic_vertex_index_map;
    };
}


