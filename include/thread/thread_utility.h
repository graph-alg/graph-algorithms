/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : ThreadUtility.h
* @brief      : common head files for thread
* @version    : 1.0
* @date       : 2019/09/02
******************************************************************************************************************/

#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <queue>
#include <thread>


#include "container/safe_queue.hpp"
#include "thread/join_thread.h"
#include "thread/smart_lock.hpp"
#include "thread/function_wrapper.hpp"

using std::atomic_uint32_t;
using std::atomic_bool;
using std::atomic_flag;
using std::condition_variable;
using std::future;
using std::mutex;
using std::queue;
using std::shared_mutex;
using std::shared_ptr;
using std::thread;

using std::for_each;
using std::mem_fn;

enum FUTURE_FLAG{
    DELETE,
    INSERT,
    IGNORE
};

