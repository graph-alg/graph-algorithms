/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : double_node.h
* @brief      : A simple bi-direct linked_node structure for storing an integer
* @version    : 1.0
* @date       : 2019/8/31
******************************************************************************************************************/

#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
      * \class A double linked linked_node for storing an integer
      */
    template <typename value_type>
    class double_node {
    public:
        /**
          * @fn double_node
          * @brief construct a new double linked linked_node with a given other_value
          * @param other_value
          */
        explicit double_node(value_type other_value): value(other_value) {

        }

        /**
         * @fn double_node
         * @brief construct a new double linked linked_node with a given double linked linked_node
         * @param node
         */
         explicit double_node(const shared_ptr<double_node>& node):
                value(node->get_value()), next(node->get_next()), prior(node->get_prior()){
        }

        /**
         * @fn at
         * @brief get the value of this linked_node
         * @return
         */
        value_type get_value() const{
            return value;
        }

        /**
         * @fn get_next
         * \brief get the next linked_node of this linked_node
         * @return
         */
        shared_ptr<double_node> get_next() const{
            return next;
        }

        /**
         * @fn get_prior
         * \brief get the prior linked_node of this linked_node
         * @return
         */
        shared_ptr<double_node> get_prior() const{
            return prior;
        }

        /**
         * @fn set_next
         * @brief set the next other_node for this other_node
         * @param other_node
         */
        void set_next(const shared_ptr<double_node>& other_node) {
            if(other_node)
            {
                next = other_node;
            }else
            {
                next.reset();
            }
        }

        /**
         * @fn set_prior
         * @brief set the prior double linked other_node for this other_node
         * @param other_node
         */
        void set_prior(const shared_ptr<double_node>& other_node) {
            if(other_node)
            {
                prior = other_node;
            }else
            {
                prior.reset();
            }
        }



    private:
        value_type value;
        shared_ptr<double_node> next;
        shared_ptr<double_node> prior;
    };
}

