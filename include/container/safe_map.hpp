/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : safe_map.hpp
* @brief      : A thread safe map
* @version    : 1.1
* @date       : 2020/10/14
******************************************************************************************************************/
#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details a thread safe map, provide STL-style functions
     * @tparam key_type
     * @tparam value_type
     */
    template <typename key_type, typename value_type>
    class safe_map
    {
    public:
        /**
         * @details get the corresponding value with the given key
         * @remarks return a copy rather than a reference
         * @param key
         * @return
         */
        value_type at(key_type key)
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.at(key);
        }

        /**
         * @details implement begin function for this map
         * @return
         */
        auto begin()
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.begin();
        }

        /**
         * @details clear data of this map
         */
        void clear()
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.clear();
        }

        /**
         * @details judge the given key in this map or not
         * @param key
         * @return
         */
        bool count(key_type key)
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.count(key);
        }

        /**
         * @details judge this map is empty or not
         */
        void empty()
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.empty();
        }

        /**
         * @details implement the end function for this map
         * @return
         */
        auto end()
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.end();
        }


        /**
         * @details erase a pair by the given key
         * @param key
         */
        void erase(key_type key)
        {
            lock_guard<mutex> lg(map_mutex);
            data_map.erase(key);
        }


        /**
         * @details insert a new pair
         * @param p
         */
        void insert(const pair<key_type,value_type>& p)
        {
            lock_guard<mutex> lg(map_mutex);
            data_map.insert(p);
        }

        /**
         * @details reset the value of the given key
         * @remarks it can ensure the security in a multi-thread environment
         * @param key
         * @param value
         */
        void set_value(key_type key, value_type value)
        {
            lock_guard<mutex> lg(map_mutex);
            data_map.at(key) = value;
        }

        /**
         * @details get the size of this map
         * @return
         */
        uint32_t size()
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.size();
        }

    private:
        mutable mutex map_mutex;
        unordered_map<key_type,value_type> data_map;
    };
}


