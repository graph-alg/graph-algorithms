/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : complicated_thread_pool.h
* @brief      : A pool of thread
* @version    : 1.0
* @date       : 2020/9/2
******************************************************************************************************************/

#pragma once
#include "thread/function_wrapper.hpp"
#include "container/safe_queue.hpp"
#include "thread/join_thread.h"
#include "thread/work_stealing_queue.h"


namespace scnu
{
    /**
     * @class complicated_thread_pool
     * @brief this pool can implement load balancing
     * @note this pool cannot be used in dependence task
     */
    class complicated_thread_pool
    {
    public:
        complicated_thread_pool();

        explicit complicated_thread_pool(uint32_t thread_number);

        ~complicated_thread_pool();

        void decrease_running_thread_count()
        {
            lock_guard<mutex> lg(running_thread_count_mutex);
            --running_thread_count;
        }

        static bool pop_task_from_local_queue(function_wrapper& type);

        bool pop_task_from_pool_queue(function_wrapper& task);

        bool pop_task_from_other_thread_queue(function_wrapper &task);

        void run_pending_task();

        template<typename function_type>
        std::future<std::result_of_t<function_type()>> submit_async_task(function_type function)
        {
            using ResultType =  std::result_of_t<function_type()>;
            std::packaged_task<ResultType()> task(function);
            std::future<ResultType> result(task.get_future());
            if(local_work_queue)
            {
                local_work_queue->push(function_wrapper(std::move(task)));
            }
            else
            {
                pool_work_queue.push(function_wrapper(std::move(task)));
            }
            return result;
        }

        template<typename function_type>
        void submit_task(function_type function)
        {
            using ResultType =  std::result_of_t<function_type()>;
            std::function<ResultType()> task(function);
            pool_work_queue.push(function_wrapper(std::move(task)));
        }

        void wait_end();

    private:
        condition_variable end_condition;

        std::atomic_bool done;

        safe_queue<function_wrapper> pool_work_queue;

        std::vector<std::unique_ptr<work_stealing_queue>> queue_vector;

        uint32_t running_thread_count{};

        mutable mutex running_thread_count_mutex;

        std::vector<std::thread> thread_vector;

        static thread_local work_stealing_queue* local_work_queue;

        static thread_local uint32_t thread_index;

        void worker_thread(uint32_t index);
    };

    thread_local work_stealing_queue* complicated_thread_pool::local_work_queue{nullptr};

    thread_local uint32_t complicated_thread_pool::thread_index{0};
}




