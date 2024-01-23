/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : core_test.h
* @brief      : A common head file for truss maintenance
* @version    : 1.0
* @date       : 2021/01/06
******************************************************************************************************************/

#pragma once

#include "container/container_compare.hpp"
#include "container/container_convert.hpp"
#include "container/container_copy.hpp"
#include "container/container_operate.hpp"
#include "truss/basic_truss_decomposition.h"
#include "truss/bin_sort_truss_decomposition.h"
#include "truss/order_truss_maintenance.h"
#include "truss/jes_truss_maintenance.h"
#include "truss/mixed_structure_truss_maintenance.h"
#include "truss/jes_order_truss_maintenance.h"
#include "truss/truss_compare.h"
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
using scnu::basic_truss_decomposition;
using scnu::bin_sort_truss_decomposition;
using scnu::container_compare;
using scnu::container_convert;
using scnu::container_copy;
using scnu::double_list;
using scnu::double_node;
using scnu::extend_list;
using scnu::extend_node;
using scnu::process_information;
using scnu::random_generator;
using scnu::simple_logger;
using scnu::simple_timer;
using scnu::truss_compare;
using scnu::order_truss_maintenance;
using scnu::jes_truss_maintenance;
using scnu::jes_order_truss_maintenance;
using scnu::mixed_structure_truss_maintenance;


