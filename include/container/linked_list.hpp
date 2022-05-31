/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : CoreMaintenance.cpp
* @brief      : a data structure for single linked list
* @version    : 1.1
* @date       : 2020/10/14
******************************************************************************************************************/

#pragma once
#include "container/linked_node.hpp"

namespace scnu
{
    /**
     * @details a single linked list for directly visiting node, provide STL style interface
     * @tparam value_type
     */
    template <typename value_type>
    class linked_list
    {
        using iterator = list_forward_iterator<shared_ptr<linked_node<value_type>>>;
    public:

        /**
         * @details get the begin iterator of this list
         * @return
         */
        iterator begin()
        {
            return iterator(head);
        }

        void clear(){
            auto p = head;
            while (p){
                auto q = p->get_next();
                p->set_next(shared_ptr<linked_node<value_type>>());
                p = q;
            }
            head = move(shared_ptr<linked_node<value_type>>());
            rear = head;
        }

        /**
         * @details get the end iterator of this list
         * @return
         */
        iterator end()
        {
            return iterator(shared_ptr<linked_node<value_type>>());
        }

        bool empty(){
            return !head ? true:false;
        }

        shared_ptr<linked_node<value_type>> emplace_back(value_type value)
        {
            auto node = make_shared<linked_node<value_type>>(move(value));
            if(!rear)
            {
                head = node;
            } else
            {
                rear->set_next(node);
            }
            rear = node;
            return node;
        }
        /**
         * @details insert new value before the rear
         * @param value
         */
        shared_ptr<linked_node<value_type>> emplace_front(value_type value)
        {
            auto node = make_shared<linked_node<value_type>>(move(value));
            if(!head)
            {
                rear = node;
            }else
            {
                node->set_next(head);
            }
            head = node;
            return node;
        }

        shared_ptr<linked_node<value_type>> erase(iterator iter){
            auto node = *iter;
            erase(node);
            return node;
        }

        void erase(const shared_ptr<linked_node<value_type>>& node){
            if(node == head){
                if(head == rear){
                    head = shared_ptr<linked_node<value_type>>();
                    rear = shared_ptr<linked_node<value_type>>();
                    head == rear;
                }else
                {
                    head = head->get_next();
                }
            }else if(node == rear)
            {
                auto prior_node = find_prior_node(node);
                rear = prior_node;
                rear->set_next(shared_ptr<linked_node<value_type>>());
            }else
            {
                auto prior_node = find_prior_node(node);
                auto next_node = node->get_next();
                prior_node->set_next(next_node);
            }
            node->set_next(shared_ptr<linked_node<value_type>>());
        }

        /**
         * @details find a node with the given value
         * @remarks it just returns the first node
         * @param value
         * @return
         */
        shared_ptr<linked_node<value_type>> find(value_type value)
        {
            auto p = head;
            while(p && p->get_value() != value)
            {
                p = p->get_next();
            }
            return p;
        }

        shared_ptr<linked_node<value_type>> find_prior_node(const shared_ptr<linked_node<value_type>>& node){
            auto p = head;
            while(p && p->get_next()!=node)
            {
                p = p->get_next();
            }
            return p;
        }
        
        /**
         * @details get the head of this list
         * @return
         */
        shared_ptr<linked_node<value_type>> get_head()
        {
            return head;
        }

        /**
         * @details get the rear of this list
         * @return
         */
        shared_ptr<linked_node<value_type>> get_rear()
        {
            return rear;
        }

        /**
         * @details connect two linked lists
         * @param other_list
         */
        void connect(const shared_ptr<linked_list<value_type>>& other_list)
        {
            if(other_list->head){
                if(!rear)
                {
                    head = move(other_list->head);
                    rear = move(other_list->rear);
                } else
                {
                    auto node = move(other_list->head);
                    rear->set_next(node);
                    rear = move(other_list->rear);
                }
                other_list->head = shared_ptr<linked_node<value_type>>();
                other_list->rear = shared_ptr<linked_node<value_type>>();
            }
        }

        shared_ptr<linked_node<value_type>> insert_after(const shared_ptr<linked_node<value_type>>& pivot_node,
                                                            value_type value)
        {
            /**
             * @details two cases: (1) pivot node does not exist
             * (2) pivot node is the rear node
             */
            if(!pivot_node||!pivot_node->get_next())
            {
                return push_back(value);
            }

            auto node = make_shared<linked_node<value_type>>(value);
            auto next_node = pivot_node->get_next();
            node->set_next(next_node);
            pivot_node->set_next(node);
            return node;
        }

        /**
         * @details insert new value before the rear
         * @param value 
         */
        shared_ptr<linked_node<value_type>> push_front(value_type value)
        {
            auto node = make_shared<linked_node<value_type>>(value);
            if(!head)
            {
                rear = node;
            }else
            {
                node->set_next(head);
            }
            head = node;
            return node;
        }

        shared_ptr<linked_node<value_type>> push_front(const shared_ptr<linked_node<value_type>> node)
        {
            if(!head)
            {
                rear = node;
            }else
            {
                node->set_next(head);
            }
            head = node;
            return node;
        }

        /**
         * @details insert new value after the rear
         * @param value
         */
        shared_ptr<linked_node<value_type>> push_back(value_type value)
        {
            auto node = make_shared<linked_node<value_type>>(value);
            if(!rear)
            {
                head = node;
            } else
            {
                rear->set_next(node);
            }
            rear = node;
            return node;
        }

        shared_ptr<linked_node<value_type>> push_back(const shared_ptr<linked_node<value_type>>& node)
        {
            if(!rear)
            {
                head = node;
            } else
            {
                rear->set_next(node);
            }
            rear = node;
            return node;
        }


        void remove(value_type value){
            auto node = find(value);
            while(node){
                erase(node);
                node = find(value);
            }
        }

    private:
        shared_ptr<linked_node<value_type>> head;
        shared_ptr<linked_node<value_type>> rear;
    };
}
