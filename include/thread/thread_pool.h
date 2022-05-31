/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : thread_pool.h
* @brief      : A pool of thread
* @version    : 1.0
* @date       : 2020/9/2
******************************************************************************************************************/

#pragma once
#include "thread/thread_utility.h"

namespace scnu
{
    /**
     * @class thread_pool
     * @brief a pool of threads
     */
    class thread_pool
    {
    public:
        thread_pool();

        explicit thread_pool(uint32_t thread_number);

        ~thread_pool();

        void barrier();

        uint32_t get_thread_id(thread::id thread_id);

        uint32_t get_thread_number();

        void wait_end();

        bool pop_task_from_local_work_queue(function_wrapper& task,uint32_t thread_id);

        void run_pending_task(uint32_t thread_id);

        template<typename function_type>
        future<std::result_of_t<function_type()>> submit_async_task(function_type function)
        {
            using ResultType =  std::result_of_t<function_type()>;
            std::packaged_task<ResultType()> task(function);
            future<ResultType> result(task.get_future());
            uint32_t index = task_id % local_work_queue_vector.size();
            ++task_id;
            local_work_queue_vector.at(index)->push(function_wrapper(move(task)));
            return result;
        }


        template<typename function_type>
        void submit_task(function_type function)
        {
            uint32_t index = task_id % local_work_queue_vector.size();
            ++task_id;
            local_work_queue_vector.at(index)->push(function_wrapper(move(function)));
        }

        template<typename function_type>
        void submit_task(function_type function, uint32_t index)
        {
            local_work_queue_vector.at(index)->push(function_wrapper(move(function)));
        }


    private:

        atomic_bool done;

        vector<shared_ptr<safe_queue<function_wrapper>>> local_work_queue_vector;

        vector<thread> thread_vector;

        unordered_map<thread::id,uint32_t> thread_id_map;

        /**
         * @brief record the count of all tasks
         * @note it is only operated by main thread
         */
        uint32_t task_id{0};

        void worker_thread(uint32_t thread_id);

    };
}




