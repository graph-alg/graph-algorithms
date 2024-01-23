/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : simple_thread_pool.h
* @brief      :
* @version    : 1.0
* @date       : 2020/10/1
******************************************************************************************************************/
#pragma once
#include "thread/thread_utility.h"
#include "thread/join_thread.h"

namespace scnu {

    class simple_thread_pool {
    public:
        simple_thread_pool();

        explicit simple_thread_pool(uint32_t thread_count);


        uint32_t get_thread_count();

        ~simple_thread_pool();

        void wait_end();

        bool pop_task_from_work_queue(function_wrapper &task);

        void run_pending_task();

        template<typename function_type>
        std::future<std::result_of_t<function_type()>> submit_async_task(function_type function) {
            using result_type = std::result_of_t<function_type()>;
            std::packaged_task<result_type()> task(function);
            std::future<result_type> result(task.get_future());
            work_queue.push(function_wrapper(std::move(task)));
            return result;
        }


        template<typename function_type>
        void submit_task(function_type function) {
            work_queue.push(function_wrapper(std::move(function)));
        }


    private:
        std::atomic_bool done;

        safe_queue<function_wrapper> work_queue;

        std::vector<std::thread> thread_vector;

        void worker_thread();

    };
}



