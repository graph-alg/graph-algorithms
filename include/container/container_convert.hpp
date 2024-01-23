/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : container_convert.hpp
* @details    : template functions for converting different containers
* @version    : 1.2
* @date       : 2020/12/9
******************************************************************************************************************/

#pragma once
#include "container/container_utility.h"
#include "container/safe_map.hpp"
#include "container/safe_set.hpp"
#include "container/safe_unordered_map.hpp"
#include "container/safe_unordered_set.hpp"

namespace scnu
{
    /**
     * @details a class for converting a container to another type of containers
     */
    class container_convert
    {
    public:
        /**
         * @details convert the given container to std map
         * @tparam key_type
         * @tparam value_type
         * @tparam source_container_type
         * @param source_container
         * @return
         */
        template<typename key_type,typename value_type,typename source_container_type>
        static shared_ptr<map<key_type,value_type>> to_map(const source_container_type& source_container)
        {
            auto destination_container = make_shared<map<key_type,value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }

        /**
        * @details convert the given container to safe map
        * @tparam key_type
        * @tparam value_type
        * @tparam source_container_type
        * @param source_container
        * @return
        */
        template<typename key_type,typename value_type,typename source_container_type>
        static shared_ptr<safe_map<key_type,value_type>> to_safe_map(const source_container_type& source_container)
        {
            auto destination_container = make_shared<safe_map<key_type,value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }

        /**
        * @details convert the given container to std set
        * @tparam value_type
        * @tparam source_container_type
        * @param source_container
        * @return
        */
        template<typename value_type,typename source_container_type>
        static shared_ptr<set<value_type>> to_set(const source_container_type& source_container)
        {
            auto destination_container = make_shared<set<value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }


        /**
         * @details convert the given container to safe set
         * @tparam value_type
         * @tparam source_container_type
         * @param source_container
         * @return
         */
        template<typename value_type,typename source_container_type>
        static shared_ptr<safe_set<value_type>> to_safe_set(const source_container_type& source_container)
        {
            auto destination_container = make_shared<safe_set<value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }

        /**
        * @details convert the given container to safe_unordered_map
        * @tparam key_type
        * @tparam value_type
        * @tparam source_container_type
        * @param source_container
        * @return
        */
        template<typename key_type,typename value_type,typename source_container_type>
        static shared_ptr<safe_unordered_map<key_type,value_type>> to_safe_unordered_map(const source_container_type& source_container)
        {
            auto destination_container = make_shared<safe_unordered_map<key_type,value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }

        /**
         * @details convert the given container to safe_unordered_set
         * @tparam value_type
         * @tparam source_container_type
         * @param source_container
         * @return
         */
        template<typename value_type,typename source_container_type>
        static shared_ptr<safe_unordered_set<value_type>> to_safe_unordered_set(const source_container_type& source_container)
        {
            auto destination_container = make_shared<safe_unordered_set<value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }

        /**
         * @details convert the given container to std unordered_map
         * @tparam key_type
         * @tparam value_type
         * @tparam source_container_type
         * @param source_container
         * @return
         */
        template<typename key_type,typename value_type,typename source_container_type>
        static shared_ptr<unordered_map<key_type,value_type>> to_unordered_map(const source_container_type& source_container)
        {
            auto destination_container = make_shared<unordered_map<key_type,value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }

        /**
         * @details convert the given container to std unordered_set
         * @tparam value_type
         * @tparam source_container_type
         * @param source_container
         * @return
         */
        template<typename value_type,typename source_container_type>
        static shared_ptr<unordered_set<value_type>> to_unordered_set(const source_container_type& source_container)
        {
            auto destination_container = make_shared<unordered_set<value_type>>();
            for(const auto& value:*source_container)
            {
                destination_container->insert(value);
            }
            return destination_container;
        }
    };
}
