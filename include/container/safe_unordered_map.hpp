/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : safe_unordered_map.hpp
* @brief      : a safe unordered map
* @version    : 1.1
* @date       : 2020/10/14
******************************************************************************************************************/
#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details a safe unordered_map, provide STL-style functions
     * @tparam key_type
     * @tparam value_type
     */
    template <typename key_type, typename value_type>
    class safe_unordered_map
    {
    public:
        /**
         * @details get the value with the given key
         * @note it returns a reference
         * @param key
         * @return
         */
        value_type& at(key_type key)
        {
            return data_map.at(key);
        }

        /**
         * @details implement begin function for this map
         * @return
         */
        auto begin()
        {
            return data_map.begin();
        }

        /**
         * @details clear the data of this map
         */
        void clear()
        {
            lock_guard<mutex> lg(map_mutex);
            data_map.clear();
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
         * @return
         */
        bool empty()
        {
            lock_guard<mutex> lg(map_mutex);
            return data_map.empty();
        }

        /**
         * @details implement end function for this map
         * @return
         */
        auto end()
        {
            return data_map.end();
        }

        /**
         * @details erase a pair from this map with the given key
         * @param key
         */
        void erase(key_type key)
        {
            lock_guard<mutex> lg(map_mutex);
            data_map.erase(key);
        }


        /**
         * @details insert a new pair into this map
         * @param p
         */
        void insert(pair<key_type,value_type> p)
        {
            lock_guard<mutex> lg(map_mutex);
            if(!data_map.count(p.first))
            {
                data_map.insert(std::move(p));
            }
        }

        void lock()
        {
            map_mutex.lock();
        }


        /**
         * @details reset the value of the given key
         * @remarks it can ensure the security in the multi-thread environment
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

        void unlock()
        {
            map_mutex.unlock();
        }

    private:
        mutable mutex map_mutex;
        unordered_map<key_type,value_type> data_map;
    };
}
