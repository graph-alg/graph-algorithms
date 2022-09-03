/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : LoggerTest.cpp
* @brief      : A test file for thread
* @version    : 1.2
* @date       : 2021/4/13
******************************************************************************************************************/

#include "container/container_compare.hpp"
#include "container/container_convert.hpp"
#include "graph/abstract_bipartite_graph.h"
#include "io/abstract_bipartite_graph_io.h"
#include "logger/simple_logger.h"
#include "random/random_generator.h"
#include "wing/basic_wing_decomposition.h"
#include "wing/wing_compare.h"
#include "wing/index_wing_decomposition.h"
#include "wing/order_wing_maintenance.h"
#include "wing/peel_wing_decomposition.h"
#include "wing/quasi_order_wing_maintenance.h"
#include "wing/quasi_wing_maintenance.h"
#include "container/extend_list.hpp"
#include "time/simple_timer.h"

using std::filesystem::path;

using scnu::LOG_RANK;

using scnu::abstract_bipartite_edge;
using scnu::abstract_bipartite_graph;
using scnu::abstract_bipartite_graph_io;
using scnu::abstract_left_vertex;
using scnu::abstract_right_vertex;
using scnu::bipartite_edge_index;
using scnu::container_compare;
using scnu::simple_logger;
using scnu::random_generator;
using scnu::basic_wing_decomposition;
using scnu::index_wing_decomposition;
using scnu::order_wing_maintenance;
using scnu::peel_wing_decomposition;
using scnu::quasi_order_wing_maintenance;
using scnu::quasi_wing_maintenance;
using scnu::extend_list;
using scnu::wing_compare;

