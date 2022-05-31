/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : ContainerUtility.hpp
* @brief      : common functions for container operate
* @version    : 1.0
* @date       : 2020/11/25
******************************************************************************************************************/

#pragma once

#include "container/container_utility.h"

namespace scnu
{
    class container_operate {
    public:
        /**
         * @details get common items of two maps
         * @tparam key_type
         * @tparam value_type
         * @param map1
         * @param map2
         * @return
         */
        template<typename key_type, typename value_type>
        static shared_ptr<map<key_type, value_type>> intersect(const shared_ptr<map<key_type, value_type>> &map1,
                                                        const shared_ptr<map<key_type, value_type>> &map2) {
            auto result_map = make_shared<map<key_type, value_type>>();
            for (const auto &[key, value]:*map1) {
                if (map2->count(key) && map2->at(key) == value) {
                    result_map->insert({key, value});
                }
            }
            return result_map;
        }

        /**
         * @details get common keys of two maps
         * @remarks items having the sub key may not equal
         * @tparam key_type
         * @tparam value_type
         * @param map1
         * @param map2
         * @return
         */
        template<typename key_type, typename value_type>
        static shared_ptr<set<key_type>> intersect_key(const shared_ptr<map<key_type, value_type>> &map1,
                                                const shared_ptr<map<key_type, value_type>> &map2) {
            auto result_set = make_shared<set<key_type>>();
            for (const auto &[key, value]:*map1) {
                if (map2->count(key)) {
                    result_set->insert(key);
                }
            }
            return result_set;
        }

        /**
         * @details get common elements of two sets
         * @tparam value_type
         * @param set1
         * @param set2
         * @return
         */
        template<typename value_type>
        static shared_ptr<set<value_type>> intersect(const shared_ptr<set<value_type>> &set1,
                                              const shared_ptr<set<value_type>> &set2) {
            auto result_set = make_shared<set<value_type>>();
            for (const auto &value:*set1) {
                if (set2->count(value)) {
                    result_set->insert(value);
                }
            }
            return result_set;
        }

        /**
         * @details get common items of two unordered_maps
         * @tparam key_type
         * @tparam value_type
         * @param map1
         * @param map2
         * @return
         */
        template<typename key_type, typename value_type>
        static shared_ptr<unordered_map<key_type, value_type>> intersect(const shared_ptr<unordered_map<key_type, value_type>> &map1,
                                                                  const shared_ptr<unordered_map<key_type, value_type>> &map2) {
            auto result_map = make_shared<unordered_map<key_type, value_type>>();
            for (const auto &[key, value]:*map1) {
                if (map2->count(key) && map2->at(key) == value) {
                    result_map->insert({key, value});
                }
            }
            return result_map;
        }

        /**
         * @details get common keys of two maps
         * @remarks items having the sub key may not equal
         * @tparam key_type
         * @tparam value_type
         * @param map1
         * @param map2
         * @return
         */
        template<typename key_type, typename value_type>
        static shared_ptr<unordered_set<key_type>> intersect_key(const shared_ptr<unordered_map<key_type, value_type>> &map1,
                                                          const shared_ptr<unordered_map<key_type, value_type>> &map2) {
            auto result_set = make_shared<unordered_set<key_type>>();
            for (const auto &[key, value]:*map1) {
                if (map2->count(key)) {
                    result_set->insert(key);
                }
            }
            return result_set;
        }

        /**
         * @details get common elements of two sets
         * @tparam value_type
         * @param set1
         * @param set2
         * @return
         */
        template<typename value_type>
        static shared_ptr<unordered_set<value_type>> intersect(const shared_ptr<unordered_set<value_type>> &set1,
                                                        const shared_ptr<unordered_set<value_type>> &set2) {
            auto result_set = make_shared<unordered_set<value_type>>();
            for (const auto &value:*set1) {
                if (set2->count(value)) {
                    result_set->insert(value);
                }
            }
            return result_set;
        }
    };
}
