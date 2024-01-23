/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : extend_list.h
* @brief      : A bi-direct node for storing the corresponding pair
* @version    : 1.2
* @date       : 2022/06/20
******************************************************************************************************************/

#pragma once
#include "container/extend_node.hpp"

namespace scnu
{


    /**
     * @details An extend double linked list for storing pairs,
     * it implements a relative order between nodes, key is regarded as a rank
     * and the value is used to distinguish nodes
     * @remarks (1) this list is not suitable for storing repeated value
     * (2) the key is automatically generate
     */
    template <typename key_type, typename value_type>
    class extend_list {
        /**
         * @brief define STL style iterator
         */
        using iterator = list_iterator<shared_ptr<extend_node<key_type,value_type>>>;
        using forward_iterator = list_forward_iterator<shared_ptr<extend_node<key_type,value_type>>>;
        using backward_iterator = list_backward_iterator<shared_ptr<extend_node<key_type,value_type>>>;
    public:
        /**
         * @remarks for complicated value_type, it should implement hash function in corresponding class
         * or store smart pointer type
         */
        extend_list():node_map(make_shared<unordered_map<value_type, shared_ptr<extend_node<key_type, value_type>>>>()) {

        }

        /**
         * @details copy a list (rear insert)
         * @param other_list
         */
        explicit extend_list(const shared_ptr<extend_list<key_type,value_type>>& other_list):extend_list(){
            auto p = other_list->get_head();
            while (p) {
                auto node = make_shared<extend_node<key_type, value_type>>(0, p->get_value());
                right_insert(node);

                p = p->get_next();
            }
        }

        /**
         * @details provide forward iterator
         * @return
         */
        iterator begin()
        {
            return forward_iterator(head);
        }


        /**
         * @details clear this list
         */
        void clear() {
            head->reset();
            rear->reset();
            node_map->clear();
        }

        /**
         * @details check this list containing the given value or not
         * @param value
         * @return
         */
        bool count_value(const value_type &value) {
            return node_map->count(value);
        }

        /**
         * @details check this list empty or not
         * @return
         */
        bool empty()
        {
            return node_map->empty();
        }

        /**
         * @details provide forward iterator
         * @return
         */
        iterator end()
        {
            return forward_iterator(shared_ptr<extend_node<key_type, value_type>>());
        }

        /**
         * @details get the head of this list
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> get_head() {
            return head;
        }

        /**
         * @details find a node with the given value
         * @remarks since the value is unique, the node is also unique
         * @param value
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> find(const value_type& value)
        {
            return node_map->count(value) ? node_map->at(value):shared_ptr<extend_node<key_type, value_type>>();
        }

        /**
         * @details find the key of a node with the given value
         * @param value
         * @return
         */
        std::optional<key_type> find_key(const value_type& value)
        {
            if(node_map->count(value))
            {
                return node_map->at(value)->get_key();
            }
            return std::nullopt;
        }

        /**
         * @details get the rear of this list
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> get_rear() {
            return rear;
        }

        /**
        * @details get the size of this list
        * @return
        */
        uint32_t size()
        {
            return node_map->size();
        }

        /**
         * @details push_back a node after the given node
         * @param node
         * @param pivot_node
         */
        void insert_after(const shared_ptr<extend_node<key_type,value_type>>& node,
                          const shared_ptr<extend_node<key_type,value_type>>& pivot_node)
        {
            if(!pivot_node)
            {
                /**
                 * @brief the pivot node does count
                 */
                head = node;
                rear = node;
                node->set_key(0);
            }else
            {
                if(pivot_node==rear)
                {
                    /**
                     * @brief the next node does not count
                     * @remarks the key is comparable and ordered
                     */
                    auto key = pivot_node->get_key() + 1;
                    node->set_key(key);

                    pivot_node->set_next(node);
                    node->set_prior(pivot_node);
                    node->set_next(nullptr);
                    rear = node;
                } else
                {
                    /**
                     * @brief ensure the new key is ordered
                     */
                    auto next_node = pivot_node->get_next();

                    pivot_node->set_next(node);
                    node->set_prior(pivot_node);
                    node->set_next(next_node);
                    next_node->set_prior(node);
                }
            }
            node_map->insert({node->get_value(), node});
        }

        /**
         * @details push_back a linked_node before the given node
         * @param node
         * @param pivot_node
         */
        void insert_before(const shared_ptr<extend_node<key_type,value_type>>& node,
                           const shared_ptr<extend_node<key_type,value_type>>& pivot_node)
        {
            if(!pivot_node)
            {
                /**
                  * @brief the pivot node does exist
                  */
                head = node;
                rear = node;
                node->set_key(0);
            } else
            {

                /**
                 * @brief the prior node does not exist
                 */
                if(pivot_node == head)
                {
                    double key = pivot_node->get_key() - 1;
                    node->set_key(key);

                    node->set_prior(nullptr);
                    node->set_next(pivot_node);

                    pivot_node->set_prior(node);
                    head = node;
                } else
                {
                    /**
                     * @brief ensure key is ordered
                     */
                    auto prior_node = pivot_node->get_prior();

                    node->set_next(pivot_node);
                    node->set_prior(prior_node);
                    prior_node->set_next(node);
                    pivot_node->set_prior(node);
                }
            }
            node_map->insert({node->get_value(), node});
        }

        /**
         * @details push_back a new node with the given value before the head of this list
         * @param value
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> left_insert(const value_type& value) {
            /**
            * @brief two cases: (1) the head of this list does not count
             * (2) the head of this list exists
            */
            auto node = make_shared<extend_node<key_type,value_type>>(0, value);
            if (!head) {
                rear = node;
            } else {
                /**
                 * @brief reset the key of this node to ensure keys are ordered
                 */
                node->set_key(head->get_key() - 1);

                node->set_next(head);
                head->set_prior(node);
            }
            head = node;
            node_map->insert({value, node});

            return node;
        }

        /**
         * @fn push_front
         * @brief push_back a new linked_node before the head of this list
         * @param node
         */
        void left_insert(const shared_ptr<extend_node<key_type,value_type>>& node) {
            /**
             * @details two cases: (1) the head of this list does not count
             * (2) the head of this list exists
             */
            if (!head) {
                node->set_key(0);
                rear = node;
            } else {
                node->set_key(head->get_key() - 1);

                node->set_next(head);
                head->set_prior(node);
            }
            head = node;
            node_map->insert({node->get_value(), node});
        }

        iterator rbegin()
        {
            return backward_iterator(rear);
        }

        iterator rend()
        {
            return backward_iterator(nullptr);
        }

        /**
          * @fn push_back
          * @brief push_back a linked_node with the given value after the rear
          * @param node
          */
        void right_insert(const shared_ptr<extend_node<key_type,value_type>>& node) {
            /**
             * @details two cases: (1) the head of this list does not count
             * (2) the head of this list exists
             */
            if (!rear) {
                node->set_key(0);
                head = node;
            } else {
                node->set_key(rear->get_key() + 1);
                rear->set_next(node);
                node->set_prior(rear);
            }
            rear = node;
            node_map->insert({node->get_value(), node});
        }

        /**
         * @fn push_back
         * @brief push_back a new linked_node with the given value
         * @param value
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> push_back(const value_type& value) {
            /**
            * @brief two cases: (1) the rear of this list does not count
            * (2) the rear of this list exists
            */
            auto node = make_shared<extend_node<key_type,value_type>>(0, value);
            if (!rear) {
                head = node;
            } else {
                node->set_key(rear->get_key() + 1);

                rear->set_next(node);
                node->set_prior(rear);
            }
            rear = node;
            node_map->insert({value, node});
            return node;
        }

        /**
         * @fn remove
         * @brief remove a linked_node with the given value
         * @param value
         * @return
         */
        shared_ptr<extend_node<key_type,value_type>> remove(const value_type& value) {
            auto node = node_map->at(value);
            if (head == rear) {
                /**
                 * @brief this list only contains one node
                 */
                head.reset();
                rear.reset();
            } else {
                if (node == head) {
                    /**
                     * @brief remove the head of this list
                     */
                    head = head->get_next();
                    head->set_prior(nullptr);
                } else if (node == rear) {
                    /**
                     * @brief remove the rear of this list
                     */
                    rear = rear->get_prior();
                    rear->set_next(nullptr);
                } else {
                    /**
                     * @brief remove the middle linked_node of this list
                     */
                    auto next_node = node->get_next();
                    auto prior_node = node->get_prior();

                    prior_node->set_next(next_node);
                    next_node->set_prior(prior_node);
                }
            }
            node_map->erase(value);

            node->set_prior(nullptr);
            node->set_next(nullptr);
            return node;
        }

        /**
         * @details set the head of this list
         * @param node
         */
        void set_head(const shared_ptr<extend_node<key_type,value_type>>& node) {
            head = node;
        }

        /**
         * @details set the rear of this list
         * @param node
         */
        void set_rear(const shared_ptr<extend_node<key_type,value_type>>& node) {
            rear = node;
        }

        void reset_left_order(const shared_ptr<extend_node<key_type,value_type>>& node){
            for(auto p = node; p ; p = p->get_prior()){
                p->set_key(p->get_next()->get_key() - 1);
            }
        }

        void reset_right_order(const shared_ptr<extend_node<key_type,value_type>>& node){
            for(auto p = node; p ; p = p->get_next()){
                p->set_key(p->get_prior()->get_key() + 1);
            }
        }

        void reset_order(){
            if(head){
                head->set_key(0);
                reset_right_order(head->get_next());
            }
        }

    private:
        shared_ptr<extend_node<key_type,value_type>> head;
        shared_ptr<extend_node<key_type,value_type>> rear;
        shared_ptr<unordered_map<value_type,shared_ptr<extend_node<key_type,value_type>>>> node_map;
    };
    
}





