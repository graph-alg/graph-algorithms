/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : truss_compare.hpp
* @brief      : A class for comparing truss number in given containers
* @version    : 1.0
* @date       : 2022/06/19
******************************************************************************************************************/

#pragma once
#include "truss/truss_utility.h"

namespace scnu
{
    /**
     * @details common functions for comparing elements in containers
     * @remarks elements in two containers should have the sub type and operator =
     */
    class truss_compare
    {
    public:
        static bool same_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container1,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container2);

        static bool sub_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container1,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container2);
    };
}

