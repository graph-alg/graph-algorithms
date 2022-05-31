
#include "thread/thread_pool.h"

namespace scnu {

    thread_pool::thread_pool() : thread_pool(std::thread::hardware_concurrency()) {

    }


    thread_pool::thread_pool(uint32_t thread_number):done(false),task_id(0)
    {
        for (uint32_t i = 0; i < thread_number; ++i) {
            local_work_queue_vector.emplace_back(make_shared<safe_queue<function_wrapper>>());
        }

        for(uint32_t i = 0; i < thread_number; ++i)
        {
            auto t = std::thread(&thread_pool::worker_thread, this,i);
            auto thread_id = t.get_id();
            thread_id_map.insert({thread_id,i});
            thread_vector.emplace_back(move(t));
        }
    }

    thread_pool::~thread_pool() {
        wait_end();
    }

    void thread_pool::barrier() {
        auto counter = make_shared<uint32_t>();
        auto counter_mutex = make_shared<mutex>();
        auto counter_cv = make_shared<condition_variable>();
        for(uint32_t i = 0;i<thread_vector.size();++i){
            submit_task([=] {
                unique_lock<mutex> lk(*counter_mutex);
                ++(*counter);
                if (*counter == thread_vector.size()) {
                    counter_cv->notify_one();
                }
            });
        }
        unique_lock<mutex> lk(*counter_mutex);
        counter_cv->wait(lk,[=]{
            return *counter == thread_vector.size();
        });

        /**
         * @brief reset task count number
         */
        task_id = 0;
    }

    uint32_t thread_pool::get_thread_id(thread::id thread_id){
        return thread_id_map.count(thread_id) ? thread_id_map.at(thread_id):0;
    }

    uint32_t thread_pool::get_thread_number()
    {
        return thread_vector.size();
    }


    bool thread_pool::pop_task_from_local_work_queue(function_wrapper &task,uint32_t thread_id) {
        auto local_work_queue = local_work_queue_vector.at(thread_id);
        return local_work_queue->try_pop(task);
    }

    void thread_pool::run_pending_task(uint32_t thread_id) {
        function_wrapper task;
        if (pop_task_from_local_work_queue(task,thread_id)) {
            task();
        } else {
            std::this_thread::yield();
        }
    }

    void thread_pool::wait_end()
    {
        while (!std::all_of(local_work_queue_vector.begin(),local_work_queue_vector.end(),[=](const auto& local_work_queue){
            return local_work_queue->empty();
        }));
        done = true;
        for_each(thread_vector.begin(),thread_vector.end(),[=](auto&t){
            if(t.joinable())
            {
                t.join();
            }
        });
    }

    void thread_pool::worker_thread(uint32_t thread_id) {
        while (!done) {
            run_pending_task(thread_id);
        }
    }
}
