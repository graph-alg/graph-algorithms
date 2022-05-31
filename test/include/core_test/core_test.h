/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : core_test.h
* @brief      : A common head file for core maintenanc
* @version    : 1.0
* @date       : 2021/01/06
******************************************************************************************************************/

#include "core/basic_core_decomposition.h"
#include "core/bin_sort_core_decomposition.h"
#include "container/container_compare.hpp"
#include "container/container_convert.hpp"
#include "container/container_copy.hpp"
#include "container/container_operate.hpp"
#include "core/core_compare.h"
#include "core/order_maintenance.h"
#include "core/jes_core_maintenance.h"
#include "core/parallel_quasi_core_maintenance.h"
#include "core/quasi_core_maintenance.h"
#include "core/order_core_maintenance.h"
#include "core/traversal_core_maintenance.h"
#include "io/abstract_graph_io.h"
#include "logger/simple_logger.h"
#include "time/simple_timer.h"
#include "random/random_generator.h"
#include "system/process_information.h"

using scnu::LOG_RANK;

using scnu::abstract_edge;
using scnu::abstract_vertex;
using scnu::abstract_graph;
using scnu::abstract_graph_io;
using scnu::basic_core_decomposition;
using scnu::bin_sort_core_decomposition;
using scnu::container_compare;
using scnu::container_convert;
using scnu::container_copy;
using scnu::core_compare;
using scnu::double_list;
using scnu::double_node;
using scnu::extend_list;
using scnu::extend_node;
using scnu::order_maintenance;
using scnu::jes_core_maintenance;
using scnu::parallel_quasi_core_maintenance;
using scnu::process_information;
using scnu::quasi_core_maintenance;
using scnu::random_generator;
using scnu::safe_set;
using scnu::safe_unordered_map;
using scnu::safe_unordered_set;
using scnu::simple_logger;
using scnu::order_core_maintenance;
using scnu::simple_timer;
using scnu::traversal_core_maintenance;


