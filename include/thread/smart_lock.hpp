/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : thread_pool.h
* @brief      : A pool of thread
* @version    : 1.0
* @date       : 2020/12/3
******************************************************************************************************************/
#pragma once
#include "thread/thread_utility.h"

namespace scnu
{
    template<typename value_type>
    class smart_lock{
    public:
        smart_lock():value_set(new unordered_set<value_type>){

        }

        void lock(const value_type& value)
        {
            while (!try_lock(value));
        }

        void lock(const std::initializer_list<value_type>& value_container)
        {
            while (!try_lock(value_container));
        }


        bool try_lock(value_type value)
        {
            lock_guard<mutex> lg(set_mutex);
            if(value_set->count(value))
            {
                return false;
            }
            value_set->insert(value);
            return true;
        }

        bool try_lock(const std::initializer_list<value_type>& value_container)
        {
            lock_guard<mutex> lg(set_mutex);
            for(const auto&value:value_container)
            {
                if(value_set->count(value))
                {
                    return false;
                }
            }
            for(const auto&value:value_container)
            {
                value_set->insert(value);
            }
            return true;
        }

        void unlock(const std::initializer_list<value_type>& value_container)
        {
            set_mutex.lock();
            for(const auto&value:value_container)
            {
                value_set->erase(value);
            }
            set_mutex.unlock();
        }

        void unlock(value_type value)
        {
            set_mutex.lock();
            value_set->erase(value);
            set_mutex.unlock();
        }


        void unlock_all()
        {
            value_set->clear();
        }

        shared_ptr<unordered_set<value_type>> get_value_set()
        {
            return value_set;
        }

    private:
        mutable mutex set_mutex;
        shared_ptr<unordered_set<value_type>> value_set;
    };
}
