/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : bin_sort_core_decomposition.h
* @brief      : A bin sort core decomposition algorithm
* @version    : 1.1
* @date       : 2020/10/16
******************************************************************************************************************/

#pragma once
#include "core/core_utility.h"

namespace scnu
{
    /**
     * @details a class of bin sort core decomposition
     */
    class bin_sort_core_decomposition
    {
    public:
        static uint32_t decompose(const std::shared_ptr<abstract_graph>& graph,
                                  const std::shared_ptr<unordered_map<uint32_t,uint32_t>>& vertex_core_map);
    };
}
