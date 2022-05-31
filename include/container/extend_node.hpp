/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : extend_node.h
* @brief      : A simple bi-direct node structure for storing a pair
* @version    : 1.0
* @date       : 2019/8/31
******************************************************************************************************************/

#pragma once

#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details A four tuples linked node structure for storing a pair.
     * By default, the key type is double type
     * @remark arbitrary two pairs have different keys and the value for any pair is unique
     */
    template <typename key_type, typename value_type>
    class extend_node{
    public:
        /**
          * @details construct a new node with given key and value
          * @param other_key
          * @param other_value
          */
        extend_node(key_type other_key, value_type other_value):
                extend_node(other_key,other_value, nullptr, nullptr){

        }

        /**
          * @details construct a new node with given tuples
          * @param other_key
          * @param other_value
          */
        extend_node(key_type other_key, value_type other_value,
                    const shared_ptr<extend_node<key_type,value_type>> &other_prior,
                    const shared_ptr<extend_node<key_type,value_type>> &other_next):
                    key(other_key),value(other_value),prior(other_prior),next(other_next){

        }

        /**
         * @details construct a new node with the give node
         * @param other_node
         */
        explicit extend_node(const shared_ptr<extend_node<key_type,value_type>> &other_node)
        : extend_node(other_node->get_key(),other_node->get_value(), other_node->get_prior(),
                      other_node->get_next()){

        }

        /**
         * @details get the key of this node
         * @return
         */
        [[nodiscard]] key_type get_key() const {
            return key;
        }


        /**
         * @details get the next node of this node
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> get_next() {
            return next;
        }

        /**
         * @details get the prior node of this node
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> get_prior() {
            return prior;
        }

        /**
         * @brief get the value of this node
         * @return
         */
        [[nodiscard]] value_type get_value() const {
            return value;
        }

        /**
         * @brief reset the key of this node to a given value
         * @param key
         */
        void set_key(key_type other_key) {
            key = other_key;
        }

        /**
         * @details set the next node to the given node
         * @remarks when the given node is nullptr, reset the next node
         * @param other_node
         */
        void set_next(const shared_ptr<extend_node<key_type,value_type>>& other_node) {
            if(other_node)
            {
                next = other_node;
            }
            else
            {
                next.reset();
            }
        }

        /**
         * @details set the prior node of this node
         * @remark when the given node is nullptr, reset the prior node
         * @param other_node
         */
        void set_prior(const shared_ptr<extend_node<key_type,value_type>>& other_node) {
            if(other_node)
            {
                prior = other_node;
            }else
            {
                prior.reset();
            }
        }

    private:
        key_type key;
        value_type value;
        shared_ptr<extend_node<key_type,value_type>> prior;
        shared_ptr<extend_node<key_type,value_type>> next;
    };
}

struct extend_node_compare
{
    bool operator()(const shared_ptr<scnu::extend_node<double,uint32_t>>& node1,
                   const shared_ptr<scnu::extend_node<double,uint32_t>>& node2) const
    {
        return node1->get_key() < node2->get_key();
    }
};
