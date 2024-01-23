
#include "thread/work_stealing_queue.h"

namespace scnu
{
    bool work_stealing_queue::empty() const
    {
        std::lock_guard<std::mutex> lg(queue_mutex);
        return data_queue.empty();
    }

    void work_stealing_queue::push(value_type data)
    {
        std::lock_guard<std::mutex> lg(queue_mutex);
        data_queue.push_front(std::move(data));
    }

    bool work_stealing_queue::try_pop(value_type &value)
    {
        std::lock_guard<std::mutex> lg(queue_mutex);
        if(data_queue.empty())
        {
            return false;
        }

        value = std::move(data_queue.front());
        data_queue.pop_front();
        return true;
    }

    bool work_stealing_queue::try_steal(value_type &value)
    {
        std::lock_guard<std::mutex> lg(queue_mutex);
        if(data_queue.empty())
        {
            return false;
        }
        value = std::move(data_queue.back());
        data_queue.pop_back();
        return true;
    }
}