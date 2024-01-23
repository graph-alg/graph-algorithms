/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : join_thread.h
* @brief      : a class for automatically joining thread
* @version    : 1.0
* @date       : 2020/9/2
******************************************************************************************************************/

#pragma once
#include "thread/thread_utility.h"

namespace scnu
{
    /**
     * @class join_thread
     * @brief a class for automatically joining thread
     */
    class join_thread
    {
    public:
        explicit join_thread(std::vector<std::thread>& other_thread_vector);

        ~join_thread();
    private:
        std::vector<std::thread>& thread_vector;
    };
}




