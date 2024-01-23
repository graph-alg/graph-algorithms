/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : thread_guard.h
* @brief      : A class for auto destruct thread
* @version    : 1.0
* @date       : 2019/09/02
******************************************************************************************************************/

#pragma once
#include "thread/thread_utility.h"

namespace scnu{
    /**
     * @class thread_guard
     * @brief a class for destruct a corresponding thread
     */
    class thread_guard {
    public:
        explicit thread_guard(std::thread& other_thread);

        ~thread_guard();

    private:
        std::thread& running_thread;
    };
}



