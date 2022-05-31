/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : io_utility.h
* @brief      : An implementation of functions for input/output
* @version    : 1.1
* @date       : 2021/09/07
******************************************************************************************************************/

#pragma once
#include <filesystem>
#include <fstream>

#include "graph/abstract_graph.h"
#include "graph/abstract_bipartite_graph.h"
#include "graph/temporal_graph.h"
#include "string/string_algorithm.h"
#include "thread/thread_pool.h"

/**
 * @details load std class
 */
using std::filesystem::path;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::shared_ptr;
using std::vector;
using std::multiset;
using std::unordered_map;


/**
 * @details load std method
 */
using std::filesystem::create_directory;
using std::filesystem::exists;
using std::make_shared;
using std::make_pair;
using std::swap;
using std::sort;
