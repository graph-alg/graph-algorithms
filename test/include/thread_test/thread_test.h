/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : LoggerTest.cpp
* @brief      : A test file for logger
* @version    : 1.0
* @date       : 2020/9/3
******************************************************************************************************************/

#pragma once

#include "thread/simple_thread_pool.h"
#include "thread/thread_pool.h"
#include "container/safe_vector.hpp"
#include "random/random_generator.h"
#include <future>
#include <iostream>
#include <mutex>

using std::unique_lock;