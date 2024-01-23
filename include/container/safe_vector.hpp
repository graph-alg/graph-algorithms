/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : safe_vector.hpp
* @brief      : A thread safe vector
* @version    : 1.0
* @date       : 2020/10/14
******************************************************************************************************************/
#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details a thread safe vector, provide STL-style functions
     * @tparam value_type
     */
    template<typename value_type>
    class safe_vector
    {
        using iterator = typename vector<value_type>::iterator;
    public:
        /**
         * @details get the element at the given index
         * @remarks return a copyable value rather than reference
         * @param index
         * @return
         */
        value_type at(uint32_t index)
        {
            lock_guard<mutex> lg(vector_mutex);
            return data_vector.at(index);
        }

        /**
         * @details implement begin function for this vector
         * @return
         */
        iterator begin()
        {
            lock_guard<mutex> lg(vector_mutex);
            return data_vector.begin();
        }

        /**
         * @details push the value into this vector
         * @remarks non-copyable way
         * @param value
         */
        void emplace_back(value_type value)
        {
            lock_guard<mutex> lg(vector_mutex);
            data_vector.emplace_back(value);
        }

        /**
         * @details judge this vector is empty or not
         * @return
         */
        bool empty()
        {
            lock_guard<mutex> lg(vector_mutex);
            return data_vector.empty();
        }

        /**
         * @details erase the element at the given position
         * @param position
         */
        void erase(iterator position)
        {
            lock_guard<mutex> lg(vector_mutex);
            return data_vector.erase(position);
        }

        /**
         * @details implement end function for this vector
         * @return
         */
        iterator end()
        {
            lock_guard<mutex> lg(vector_mutex);
            return data_vector.end();
        }

        /**
         * @details clear the data of this vector
         */
        void clear()
        {
            lock_guard<mutex> lg(vector_mutex);
            data_vector.clear();
        }

        /**
         * @details push the value into this vector
         * @param value
         */
        void push_back(value_type value)
        {
            lock_guard<mutex> lg(vector_mutex);
            data_vector.push_back(value);
        }

        /**
         * @details set the new value to the given index
         * @remarks it is suitable for revising elements
         * @param index
         * @param value
         */
        void set_value(uint32_t index, value_type value)
        {
            lock_guard<mutex> lg(vector_mutex);
            data_vector.at(index) = value;
        }

        /**
         * @details get the size of this vector
         * @return
         */
        uint32_t size()
        {
            lock_guard<mutex> lg(vector_mutex);
            return data_vector.size();
        }

    private:
        mutable mutex vector_mutex;
        vector<value_type> data_vector;
    };
}