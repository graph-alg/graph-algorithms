/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : list_backward_iterator.hpp
* @brief      :
* @version    : 1.0
* @date       : 2020/10/12
******************************************************************************************************************/
#pragma once
#include "container/list_iterator.hpp"

namespace scnu
{
    /**
     * @details a derived class for implement bi-direction iterators
     * @tparam node_type
     */
    template <typename node_type>
    class list_backward_iterator: public  list_iterator<node_type>{
    public:
        /**
         * @details initialize an iterator
         * @param position
         */
        explicit list_backward_iterator(const node_type& position):
                list_iterator<node_type>(position)
        {

        }

        /**
         * @details define the incremental method of this iterator
         * @return
         */
        list_iterator<node_type>& operator++()
        {
            list_iterator<node_type>::iter = list_iterator<node_type>::iter->get_prior();
            return *this;
        }
    };


    /**
     * @class a class for implement range-based loop
     * @tparam node_type
     */
    template<typename node_type>
    class list_backward_iterator_range
    {
    public:
        /**
         * @details initialize the range
         * @param other_from
         * @param other_to
         */
        list_backward_iterator_range(const node_type& other_from,
                                    const node_type& other_to):
                from (other_from), to(other_to)
        {

        }

        /**
         * @details get the begin function
         * @return
         */
        list_backward_iterator<node_type> begin()
        {
            return from;
        }

        /**
         * @details get the end function
         * @return
         */
        list_backward_iterator<node_type> end()
        {
            return to;
        }
    private:
        list_backward_iterator<node_type> from;
        list_backward_iterator<node_type> to;
    };
}


