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

        template<class container_type>
        auto split_task(const container_type& container){
            auto thread_number = get_thread_number();
            auto location_vector = make_shared<vector<shared_ptr<decltype(container->begin())>>>(thread_number + 1,
                                                                                                 shared_ptr<decltype(container->begin())>());
            location_vector->at(0) = make_shared<decltype(container->begin())>(container->begin());
            location_vector->at(thread_number) = make_shared<decltype(container->end())>(container->end());
            if(container->empty()){
                for (uint32_t i = 1; i < thread_number; ++i) {
                    location_vector->at(i) = location_vector->at(thread_number);
                }
            }else{
                auto task_count = container->size() / thread_number + 1;
                for (uint32_t i = thread_number - 1; i > 0;--i) {
                    submit_task([=] {
                        auto distance = i * task_count;
                        if (distance < container->size()) {
                            location_vector->at(i) = make_shared<decltype(container->begin())>(container->begin());
                            std::advance(*location_vector->at(i), distance);
                        } else {
                            location_vector->at(i) = location_vector->at(thread_number);
                        }
                    });
                }
                barrier();
            }
            return location_vector;
        }

        template<typename function_type>
        future<std::result_of_t<function_type()>> submit_async_task(function_type function)
        {
            using ResultType =  std::result_of_t<function_type()>;
            std::packaged_task<ResultType()> task(function);
            future<ResultType> result(task.get_future());
            uint32_t index = task_id % local_work_queue_vector.size();
            ++task_id;
            local_work_queue_vector.at(index)->push(function_wrapper(std::move(task)));
            return result;
        }


        template<typename function_type>
        void submit_task(function_type function)
        {
            uint32_t index = task_id % local_work_queue_vector.size();
            ++task_id;
            local_work_queue_vector.at(index)->push(function_wrapper(std::move(function)));
        }

        template<typename function_type>
        void submit_task(function_type function, uint32_t index)
        {
            local_work_queue_vector.at(index)->push(function_wrapper(std::move(function)));
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




