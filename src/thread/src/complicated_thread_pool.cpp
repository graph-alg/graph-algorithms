
#include "thread/complicated_thread_pool.h"

namespace scnu {
    complicated_thread_pool::complicated_thread_pool() : complicated_thread_pool(std::thread::hardware_concurrency()) {

    }

    complicated_thread_pool::complicated_thread_pool(uint32_t thread_number) : done(false),running_thread_count(thread_number)
    {
        try {
            for (uint32_t i = 0; i < thread_number; ++i) {
                queue_vector.push_back(std::make_unique<work_stealing_queue>());
            }
            for(uint32_t i = 0; i < thread_number; ++i)
            {
                thread_vector.emplace_back(std::thread(&complicated_thread_pool::worker_thread, this, i));
            }
        }
        catch (...) {
            done = true;
            throw;
        }
    }

    complicated_thread_pool::~complicated_thread_pool() {
        wait_end();
    }

    bool complicated_thread_pool::pop_task_from_local_queue(function_wrapper &type) {
        return local_work_queue && local_work_queue->try_pop(type);
    }

    bool complicated_thread_pool::pop_task_from_pool_queue(function_wrapper &task) {
        return pool_work_queue.try_pop(task);
    }

    bool complicated_thread_pool::pop_task_from_other_thread_queue(function_wrapper &task) {
        for (unsigned int i = 0; i < queue_vector.size(); ++i) {
            unsigned const index = (thread_index + i + 1) % queue_vector.size();
            if (queue_vector.at(index)->try_steal(task)) {
                return true;
            }
        }
        return false;
    }

    void complicated_thread_pool::run_pending_task() {
        function_wrapper task;
        if (pop_task_from_local_queue(task) || pop_task_from_pool_queue(task)
            || pop_task_from_other_thread_queue(task)) {
            task();
        } else {
            decrease_running_thread_count();
            end_condition.notify_one();
            std::this_thread::yield();
        }
    }

    void complicated_thread_pool::wait_end()
    {
        unique_lock<mutex> ul(running_thread_count_mutex);
        end_condition.wait(ul,[this]{return running_thread_count==0;});
        ul.unlock();
        for_each(thread_vector.begin(),thread_vector.end(),mem_fn(&thread::join));
        done = true;
    }

    void complicated_thread_pool::worker_thread(uint32_t index) {
        thread_index = index;
        local_work_queue = queue_vector.at(thread_index).get();
        while (!done) {
            run_pending_task();
        }
    }
}