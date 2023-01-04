/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : safe_queue.h
* @brief      : A thead safe queue
* @version    : 1.1
* @date       : 2020/10/14
******************************************************************************************************************/

#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details a thread safe queue, provide STl-style functions
     * @tparam value_type
     */
    template<typename value_type>
    class safe_queue {
    public:
        /**
         * @details initialize an empty safe queue
         */
        safe_queue() = default;


        /**
         * @details initialize a thread safe queue by the given queue
         * @param other_safe_queue
         */
        safe_queue(const safe_queue& other_safe_queue)
        {
            lock_guard<mutex> lg(queue_mutex);
            data_queue = other_safe_queue.data_queue;
        }
        
        safe_queue& operator=(const safe_queue&)= delete;

        /**
         * @details judge this queue is empty or not
         * @return
         */
        bool empty() const
        {
            lock_guard<mutex> lg(queue_mutex);
            return data_queue.empty();
        }

        /**
         * @details push a new value into this queue
         * @param new_value
         */
        void push(value_type new_value)
        {
            lock_guard<mutex> lg(queue_mutex);
            data_queue.push(move(new_value));
            mDataCondition.notify_one();
        }

        /**
         * @details pop a value from this queue
         * @param value
         * @return
         */
        bool try_pop(value_type& value)
        {
            lock_guard<mutex> lg(queue_mutex);
            if(data_queue.empty())
            {
                return false;
            }
            value = std::move(data_queue.front());
            data_queue.pop();
            return true;
        }

        /**
         * @details  pop a value from this queue
         * @return
         */
        shared_ptr<value_type> try_pop()
        {
            lock_guard<mutex> lg(queue_mutex);
            if(data_queue.empty())
            {
                return shared_ptr<value_type>();
            }
            auto result = make_shared<value_type>(data_queue.front());
            data_queue.pop();
            return result;
        }

        /**
         * @details wait and pop a value from this queue
         * @param value
         */
        void wait_and_pop(value_type& value)
        {
            unique_lock<mutex> ul(queue_mutex);
            mDataCondition.wait(ul, [this]{ return !data_queue.empty();});
            value = data_queue.front();
            data_queue.pop();
            ul.unlock();
        }

        /**
         * @details wait and pop a value from this queue
         * @param value
         */
        shared_ptr<value_type> wait_and_pop()
        {
            unique_lock<mutex> ul(queue_mutex);
            mDataCondition.wait(ul, [this]{ return !data_queue.empty();});
            auto result = make_shared<value_type>(data_queue.front());
            data_queue.pop();
            ul.unlock();
            return result;
        }

    private:
        mutable mutex queue_mutex;
        queue<value_type> data_queue;
        condition_variable mDataCondition;
    };
};



