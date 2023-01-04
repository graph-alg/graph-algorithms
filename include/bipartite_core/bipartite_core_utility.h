/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bipartite_core_utility.h
* @details    : common head files for bipartite core
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/

#pragma once
#include "container/container_copy.hpp"
#include "container/double_list.hpp"
#include "container/extend_list.hpp"
#include "container/linked_list.hpp"
#include "container/safe_map.hpp"
#include "container/safe_unordered_map.hpp"
#include "container/safe_unordered_set.hpp"
#include "graph/abstract_bipartite_graph.h"
#include "thread/thread_pool.h"
#include "time/simple_timer.h"

using std::back_inserter;
using std::copy;
using std::copy_if;
using std::max;
using std::min;
using std::swap;

using scnu::abstract_bipartite_graph;
using scnu::container_copy;
using scnu::thread_pool;



