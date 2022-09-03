/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : multiple_core_test.h
* @brief      : A common head file for temporal core decomposition and maintenance
* @version    : 1.0
* @date       : 2021/09/18
******************************************************************************************************************/

#pragma once
#include "core/basic_core_decomposition.h"
#include "core/bin_sort_core_decomposition.h"
#include "container/container_compare.hpp"
#include "container/container_convert.hpp"
#include "container/container_copy.hpp"
#include "container/container_operate.hpp"
#include "io/temporal_graph_io.h"
#include "multiple_core/basic_multiple_core_decomposition.h"
#include "multiple_core/branch_multiple_core_decomposition.h"
#include "multiple_core/hierarchy_multiple_core_decomposition.h"
#include "multiple_core/quasi_multiple_core_maintenance.h"
#include "multiple_core/multiple_core.h"
#include "multiple_core/multiple_core_compare.h"
#include "multiple_core/multiple_core_query.hpp"
#include "logger/simple_logger.h"
#include "time/simple_timer.h"
#include "random/random_generator.h"


using scnu::LOG_RANK;

using scnu::basic_multiple_core_decomposition;
using scnu::container_compare;
using scnu::container_convert;
using scnu::container_copy;
using scnu::branch_multiple_core_decomposition;
using scnu::hierarchy_multiple_core_decomposition;
using scnu::random_generator;
using scnu::simple_logger;
using scnu::simple_timer;
using scnu::temporal_edge;
using scnu::multiple_core;
using scnu::multiple_core_compare;
using scnu::temporal_graph;
using scnu::temporal_graph_io;
using scnu::temporal_vertex;
using scnu::multiple_core_pair_map_index;
using scnu::multiple_core_pair_set_index;
using scnu::multiple_core_number_set_index;
using scnu::multiple_core_query;

