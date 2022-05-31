
#include "thread/join_thread.h"

namespace scnu
{
    join_thread::join_thread(std::vector<std::thread>& other_thread_vector):
            thread_vector(other_thread_vector)
    {

    }

    join_thread::~join_thread()
    {
        for(auto & t : thread_vector)
        {
            if(t.joinable())
            {
                t.join();
            }
        }
    }
}
