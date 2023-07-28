/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : truss_utility.h
* @brief      : common header files for k-truss
* @version    : 1.0
* @date       : 2022/6/17
******************************************************************************************************************/

#pragma once

#include <iterator>

#include "container/container_copy.hpp"
#include "container/double_list.hpp"
#include "container/extend_list.hpp"
#include "container/linked_list.hpp"
#include "graph/abstract_graph.h"
#include "thread/thread_pool.h"
#include "time/simple_timer.h"

using std::back_inserter;
using std::copy;
using std::copy_if;
using std::min;
using std::swap;


using scnu::abstract_edge;
using scnu::abstract_graph;
using scnu::abstract_vertex;
using scnu::container_copy;
using scnu::extend_node;
using scnu::extend_list;
using scnu::thread_pool;
using scnu::simple_timer;





