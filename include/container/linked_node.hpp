/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : linked_node.hpp
* @brief      : a linked_node structure with a single pointer
* @version    : 1.0
* @date       : 2020/9/6
******************************************************************************************************************/

#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details a linked node for storing value
     * @tparam value_type
     */
    template<typename value_type>
    class linked_node
    {
    public:
        explicit linked_node(value_type other_value): value(other_value)
        {

        }

        explicit linked_node(const shared_ptr<linked_node<value_type>>& other_node):
        value(other_node->value),next(other_node->next)
        {

        }

        ~linked_node()
        {
            next = shared_ptr<linked_node<value_type>>();
        }

        shared_ptr<linked_node<value_type>> get_next()
        {
            return next;
        }

        value_type get_value()
        {
            return value;
        }

        void set_next(const shared_ptr<linked_node>& next_node)
        {
            next = next_node;
        }

    private:
        value_type value;
        shared_ptr<linked_node> next;
    };
}


