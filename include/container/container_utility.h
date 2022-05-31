/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : container_utility.hpp
* @brief      : common header files for container
* @version    : 1.0
* @date       : 2020/11/24
******************************************************************************************************************/

#pragma once

#include <algorithm>
#include <condition_variable>
#include <iterator>
#include <mutex>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::condition_variable;
using std::lock_guard;
using std::map;
using std::mutex;
using std::nullopt;
using std::optional;
using std::pair;
using std::queue;
using std::set;
using std::shared_ptr;
using std::unordered_map;
using std::unordered_set;
using std::unique_lock;
using std::vector;

using std::copy;
using std::make_optional;
using std::make_shared;
using std::move;
using std::inserter;
using std::transform;

