
#include "thread/thread_guard.h"

namespace scnu
{
    /**
     * @fn thread_guard
     * @brief construct a thread guard
     * @param other_thread
     */
    thread_guard::thread_guard(std::thread& other_thread): running_thread(other_thread){

    }

    /**
     * @fn thread_guardd
     * @brief destruct a thread guard
     */
    thread_guard::~thread_guard()
    {
        if(running_thread.joinable())
        {
            running_thread.join();
        }
    }
}
