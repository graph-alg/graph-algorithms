/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : work_stealing_queue.h
* @brief      :
* @version    : 1.0
* @date       : 2020/9/4
******************************************************************************************************************/

#pragma once
#include "thread/function_wrapper.hpp"

namespace scnu
{
    typedef function_wrapper value_type;
    class work_stealing_queue
    {
    public:
        work_stealing_queue()=default;

        work_stealing_queue(const work_stealing_queue&) = delete;

        work_stealing_queue& operator=(const work_stealing_queue&) = delete;

        bool empty() const;

        void push(value_type value);

        bool try_pop(value_type& value);

        bool try_steal(value_type& value);

    private:
        std::deque<value_type> data_queue;

        mutable std::mutex queue_mutex;
    };
}




