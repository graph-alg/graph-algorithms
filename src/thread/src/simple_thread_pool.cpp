
#include "thread/simple_thread_pool.h"

namespace scnu
{
    simple_thread_pool::simple_thread_pool() : simple_thread_pool(std::thread::hardware_concurrency()) {

    }

    simple_thread_pool::simple_thread_pool(uint32_t thread_count) : done(false)
    {
        try {

            for(uint32_t i = 0; i < thread_count; ++i)
            {
                auto t = std::thread(&simple_thread_pool::worker_thread, this);
                thread_vector.emplace_back(move(t));
            }
        }
        catch (...) {
            done = true;
            throw;
        }
    }


    uint32_t simple_thread_pool::get_thread_count() {
        return thread_vector.size();
    }

    simple_thread_pool::~simple_thread_pool() {
        wait_end();
    }


    void simple_thread_pool::wait_end()
    {
        while (!work_queue.empty());
        done = true;
        for_each(thread_vector.begin(), thread_vector.end(),
                 [](auto &t) {
                     if (t.joinable()) {
                         t.join();
                     }
                 });
    }

    bool simple_thread_pool::pop_task_from_work_queue(function_wrapper &task) {
        return work_queue.try_pop(task);
    }

    void simple_thread_pool::run_pending_task() {
        function_wrapper task;
        if (pop_task_from_work_queue(task)) {
            task();
        } else {
            std::this_thread::yield();
        }
    }

    void simple_thread_pool::worker_thread() {
        while (!done) {
            run_pending_task();
        }
    }
}
