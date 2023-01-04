/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_core_degree_index.h
* @details    : auxiliary data structure for bipartite core
* @version    : 1.0
* @date       : 2022/10/05
******************************************************************************************************************/

#pragma once
#include "bipartite_core/bipartite_core_utility.h"

namespace scnu{
    class bipartite_core_degree_index {
    public:
        bipartite_core_degree_index();

        bipartite_core_degree_index(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_left_map,
                                    const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_right_map);

        explicit bipartite_core_degree_index(const shared_ptr<bipartite_core_degree_index>& other_index);

        bool count_left_vertex(uint32_t l);

        bool count_right_vertex(uint32_t r);

        bool empty();

        uint32_t get_left_degree(uint32_t l);

        uint32_t get_right_degree(uint32_t r);


        shared_ptr<unordered_map<uint32_t, uint32_t>> get_left_map();

        shared_ptr<unordered_map<uint32_t, uint32_t>> get_right_map();

        void insert_left_degree(uint32_t l, uint32_t degree);

        void insert_right_degree(uint32_t r, uint32_t degree);

        void insert_left_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_left_map);

        void insert_right_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_right_map);

        void update_left_degree(uint32_t l, uint32_t degree);

        void update_right_degree(uint32_t r, uint32_t degree);

        void update_left_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_left_map);

        void update_right_degree_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& other_right_map);

    private:
        shared_ptr<unordered_map<uint32_t, uint32_t>> left_map;
        shared_ptr<unordered_map<uint32_t, uint32_t>> right_map;
    };
};


