/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : list_iterator.h
* @brief      :
* @version    : 1.0
* @date       : 2020/10/12
******************************************************************************************************************/

#pragma once
#include "container/container_utility.h"
namespace scnu {
    /**
     * @details a basic list iterator for implement STL-style iterator
     * @tparam node_type
     */
    template<typename node_type>
    class list_iterator {
    public:
        /**
         * @details initialize the iterator
         * @param position
         */
        explicit list_iterator(const node_type &position) :
                iter(position) {

        }

        /**
         * @details virtual destruct function
         */
        virtual ~list_iterator() = default;

        /**
         * @details overload * operator
         * @return
         */
        auto operator*() {
            return iter;
        }

        /**
         * @details a virtual function for overloading ++ operator
         * @return
         */
        virtual list_iterator& operator++() = 0;


        /**
         * @details overload != operator
         * @param other
         * @return
         */
        bool operator!=(const list_iterator &other) const {
            return iter != other.iter;
        }

        /**
         * details over == operator
         * @param other
         * @return
         */
        bool operator==(const list_iterator &other) const {
            return iter == other.iter;
        }

    protected:
        node_type iter;
    };

    template<typename node_type>
    class list_iterator_range
    {
    public:
        list_iterator_range(const node_type& other_from,
                                    const node_type& other_to):
                from (other_from), to(other_to)
        {

        }

        list_iterator<node_type> begin()
        {
            return from;
        }

        list_iterator<node_type> end()
        {
            return to;
        }
    private:
        list_iterator<node_type> from;
        list_iterator<node_type> to;
    };
}

/**
 * @details concrete this iterator for implement STL-style iterator
 */
namespace std {
    template <typename node_type>
    struct iterator_traits<scnu::list_iterator<node_type>> {
        using iterator_category = std::forward_iterator_tag;
        using value_type = node_type;
    };
}



