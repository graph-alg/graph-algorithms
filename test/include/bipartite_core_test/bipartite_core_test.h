/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : CoreUtility.h
* @brief      : common header files for bipartite core
* @version    : 1.0
* @date       : 2021/02/01
******************************************************************************************************************/


#include "bipartite_core/bipartite_core_compare.h"
#include "bipartite_core/basic_bipartite_core_decomposition.h"
#include "bipartite_core/bipartite_core_order_index.h"
#include "bipartite_core/bipartite_core_rem_degree_index.h"
#include "bipartite_core/bipartite_core_degree_index.h"
#include "bipartite_core/branch_bipartite_core_decomposition.h"
#include "bipartite_core/branch_bipartite_core_maintenance.h"
#include "bipartite_core/share_bipartite_core_decomposition.h"
#include "bipartite_core/edge_bipartite_core_maintenance.h"
#include "bipartite_core/quasi_bipartite_core_maintenance.h"
#include "container/container_compare.hpp"
#include "container/container_convert.hpp"
#include "container/container_copy.hpp"
#include "container/container_operate.hpp"
#include "io/abstract_bipartite_graph_io.h"
#include "logger/simple_logger.h"
#include "random/random_generator.h"
#include "system/process_information.h"
#include "time/simple_timer.h"

using scnu::LOG_RANK;
using scnu::abstract_bipartite_edge;
using scnu::abstract_left_vertex;
using scnu::abstract_right_vertex;
using scnu::abstract_bipartite_graph;
using scnu::abstract_bipartite_graph_io;
using scnu::basic_bipartite_core_decomposition;
using scnu::bipartite_core;
using scnu::bipartite_core_order_index;
using scnu::bipartite_core_rem_degree_index;
using scnu::bipartite_core_degree_index;
using scnu::bipartite_core_compare;
using scnu::bipartite_core_degree_index;
using scnu::branch_bipartite_core_decomposition;
using scnu::branch_bipartite_core_maintenance;
using scnu::container_compare;
using scnu::container_convert;
using scnu::container_copy;
using scnu::double_list;
using scnu::double_node;
using scnu::share_bipartite_core_decomposition;
using scnu::edge_bipartite_core_maintenance;
using scnu::extend_list;
using scnu::extend_node;
using scnu::bipartite_core_left_store_index;
using scnu::bipartite_core_right_store_index;
using scnu::process_information;
using scnu::quasi_bipartite_core_maintenance;
using scnu::random_generator;
using scnu::safe_set;
using scnu::safe_unordered_map;
using scnu::safe_unordered_set;
using scnu::simple_logger;
using scnu::simple_timer;


