/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : double_list.h
* @details    : A bi-direction linked list without head for storing integer nodes.
* @version    : 1.0
* @date       : 2020/8/31
******************************************************************************************************************/

#pragma once
#include "container/double_node.hpp"
#include "container/list_forward_iterator.hpp"
#include "container/list_backward_iterator.hpp"

namespace scnu
{
    /**
     * @Details compared with STL list, it provides some extra operations
     * @remarks this is not a stl-style container, so there are no `begin` and `end` function
     */
    template<typename value_type>
    class double_list {
        using backward_iterator = list_backward_iterator<shared_ptr<double_node<value_type>>>;
        using forward_iterator = list_forward_iterator<shared_ptr<double_node<value_type>>>;
        using iterator = list_iterator<shared_ptr<double_node<value_type>>>;
    public:
        double_list():head(shared_ptr<double_node<value_type>>()),rear(shared_ptr<double_node<value_type>>()){

        }

        forward_iterator begin()
        {
            return forward_iterator(head);
        }
        /**
         * @details clear this list and reset head and rear to null pointer
         * @remarks this function does not release memory
         */
        void clear() {
            head.reset();
            rear.reset();
        }

        /**
         * @details judge this list is empty or not
         * @return
         */
        bool empty() {
            if(head)
            {
                return false;
            }
            return true;
        }

        forward_iterator end()
        {
            return forward_iterator(shared_ptr<double_node<value_type>>());
        }

        /**
         * @details find a node with the given value
         * @param value
         * @return
         */
        shared_ptr<double_node<value_type>> find(value_type value)
        {
            auto node = head;
            while (node)
            {
                if(node->get_value() == value)
                {
                    break;
                }
                node = node->get_next();
            }
            return node;
        }

        /**
         * @details get the head linked_node of this list
         * @return
         */
        shared_ptr<double_node<value_type>> get_head() {
            return head;
        }

        /**
         * @details get the rea linked_node of this list
         * @return
         */
        shared_ptr<double_node<value_type>> get_rear() {
            return rear;
        }

        /**
         * @details push_back a new node with the given value before the head
         * @param value
         * @return
         */
        shared_ptr<double_node<value_type>> left_insert(value_type value) {
            auto node = make_shared<double_node<value_type>>(value);
            if (!head) {
                rear = node;
            } else {
                node->set_next(head);
                head->set_prior(node);
            }
            head = node;
            return node;
        }

        /**
         * @details push_back a new double node before the head
         * @param node
         */
        void left_insert(const shared_ptr<double_node<value_type>>& node) {
            /**
            * @remark ensure nodes do not point other double node
            */
            node->set_prior(nullptr);
            node->set_next(nullptr);

            if (!head) {
                rear = node;
            } else {
                node->set_next(head);
                head->set_prior(node);
            }
            head = node;
        }

        backward_iterator rbegin()
        {
            return backward_iterator(rear);
        }

        backward_iterator rend()
        {
            return backward_iterator(shared_ptr<double_node<value_type>>());
        }

        /**
         * @details push_back a new double node with the given value after the rear
         * @param value
         * @return
         */
        shared_ptr<double_node<value_type>> right_insert(uint32_t value) {
            auto node = make_shared<double_node<value_type>>(value);
            if (!rear) {
                head = node;
            } else {
                rear->set_next(node);
                node->set_prior(rear);
            }
            rear = node;
            return node;
        }

        /**
         * @details push_back a new double node after the rear
         * @param node
         */
        void right_insert(const shared_ptr<double_node<value_type>>& node) {
            /**
             * @remark ensure nodes do not point to other double node
             */
            node->set_prior(nullptr);
            node->set_next(nullptr);

            if (!rear) {
                head = node;
            } else {
                rear->set_next(node);
                node->set_prior(rear);
            }
            rear = node;
        }

        /**
         * @fn push_back
         * @details push_back a series of nodes after the rear
         * @param start
         * @param end
         */
        void right_insert(const shared_ptr<double_node<value_type>>& start,
                          const shared_ptr<double_node<value_type>>& end) {
            /**
             * @remark ensure nodes do not point other double node
             */
             start->set_prior(nullptr);
             end->set_next(nullptr);

            if (!rear) {
                head = start;
            } else {
                rear->set_next(start);
                start->set_prior(rear);
            }
            rear = end;
        }

        /**
         * @details reset the rear to a given double node
         * @param node
         */
        void set_rear(const shared_ptr<double_node<value_type>>& node) {
            rear = node;
        }

        /**
         * @details reset the head to a given linked_node
         * @param node
         */
        void set_head(const shared_ptr<double_node<value_type>>& node) {
            head = node;
        }

        /**
         * @details remove the give linked_node from this list,
         * there are three cases: (1) head (2) rear (3) middle
         * @param node
         */
        void remove(const shared_ptr<double_node<value_type>>& node) {
            if (head == rear) {
                /**
                 * @brief the list only contains one node
                 */
                this->head.reset();
                this->rear.reset();
            } else {
                if (node == head) {
                    /**
                     * @brief delete the head node
                     */
                    head = head->get_next();
                    head->set_prior(nullptr);
                } else if (node == rear) {
                    /**
                     * @brief delete the rear node
                     */
                    rear = rear->get_prior();
                    rear->set_next(nullptr);
                } else {
                    /**
                     * @brief delete the middle linked_node
                     */
                    auto next_node = node->get_next();
                    auto prior_node = node->get_prior();
                    prior_node->set_next(next_node);
                    next_node->set_prior(prior_node);
                }
            }
        }

         /**
          * @details remove the given value from the list
          * @remarks this list only removes the value one time
          * @param value
          * @return
          */
        bool remove(value_type value) {
            auto node = find(value);
            /**
             * @details cannot find a node with the given value
             */
            if(!node)
            {
                return false;
            }
            remove(node);
            return true;
        }

    private:
        shared_ptr<double_node<value_type>> head;
        shared_ptr<double_node<value_type>> rear;
    };
}


