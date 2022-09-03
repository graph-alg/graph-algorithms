/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : graph_utility.h
* @brief      : A common head files for graph structure
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/

#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

using std::hash;
using std::map;
using std::list;
using std::nullopt;
using std::optional;
using std::set;
using std::stack;
using std::string;
using std::stringstream;
using std::shared_ptr;
using std::tuple;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using std::copy;
using std::inserter;
using std::make_pair;
using std::make_optional;
using std::make_shared;
using std::make_tuple;
using std::make_unique;
using std::swap;


/**
 * @details compare the order of two pairs
 */
struct compare_pair {
    bool operator()(const std::pair<uint64_t, uint64_t> &p1, const std::pair<uint64_t, uint64_t> &p2) const {
        return p1.second < p2.second;
    }
};

/**
 * @details the equal struct for a pair of unsigned integer value
 */
struct equal_pair {
    bool operator()(const std::pair<uint64_t, uint64_t> &p1, const std::pair<uint64_t, uint64_t> &p2) const {
        return p1.first == p2.first && p1.second == p2.second;
    }
};

/**
 * @details the hash struct for a pair of unsigned integer value
 */
struct hash_pair {
    size_t operator()(const std::pair<uint64_t , uint64_t> &p) const {
        return hash<uint64_t>()( (p.first << 2) - p.first + p.second);
    }
};



