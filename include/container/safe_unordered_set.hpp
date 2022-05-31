/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : safe_unordered_set.hpp
* @brief      : A thread safe unordered set
* @version    : 1.0
* @date       : 2020/9/16
******************************************************************************************************************/
#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details a thread safe unordered set, provide STL-style functions
     * @tparam value_type
     */
    template<typename value_type>
    class safe_unordered_set
    {
    public:
        /**
         * @details implement begin function for this set
         * @return
         */
        auto begin()
        {
            return data_set.begin();
        }

        /**
         * @details clear the data of this set
         */
        void clear()
        {
            lock_guard<mutex> lg(set_mutex);
            data_set.clear();
        }

        /**
         * @details judge the given value in the set or not
         * @param value
         * @return
         */
        bool count(value_type value)
        {
            lock_guard<mutex> lg(set_mutex);
            return data_set.count(value);
        }

        /**
         * @details judge this set is empty or not
         * @return
         */
        bool empty()
        {
            lock_guard<mutex> lg(set_mutex);
            return data_set.empty();
        }

        /**
         * @details implement end function for this set
         * @return
         */
        auto end()
        {
            return data_set.end();
        }

        /**
         * @details erase the value from this set
         * @param value
         */
        void erase(value_type value)
        {
            lock_guard<mutex> lg(set_mutex);
            data_set.erase(value);
        }

        /**
         * @details insert a new value into this set
         * @param value
         */
        void insert(value_type value)
        {
            lock_guard<mutex> lg(set_mutex);
            data_set.insert(value);
        }

        /**
         * @details get the size of this set
         * @return
         */
        uint32_t size()
        {
            lock_guard<mutex> lg(set_mutex);
            return data_set.size();
        }

    private:
        mutable mutex set_mutex;
        unordered_set<value_type> data_set;
    };
}
