/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : list_forward_iterator.hpp
* @brief      :
* @version    : 1.0.0.0
* @date       : 2020/10/12
******************************************************************************************************************/
#pragma once
#include "container/list_iterator.hpp"

namespace scnu
{
    template <typename node_type>
    class list_forward_iterator: public  list_iterator<node_type>{
    public:
        explicit list_forward_iterator(const node_type& position):
                list_iterator<node_type>(position)
        {

        }

        list_forward_iterator& operator++()
        {
            list_iterator<node_type>::iter = list_iterator<node_type>::iter->get_next();
            return *this;
        }
    };
}

