/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : HashMapCompare.hpp
* @brief      : A template class for copy elements between containers
* @version    : 1.0
* @date       : 2020/01/17
******************************************************************************************************************/

#pragma once
#include "container/container_utility.h"

namespace scnu{
    class container_copy{
    public:
        template<typename key_type,typename value_type, typename source_container_type>
        static shared_ptr<map<key_type,value_type>> to_map(const source_container_type& source_container){
            auto destination_container = make_shared<map<key_type,value_type>>();
            destination_container->reserve(source_container->size());
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container, destination_container->end()));
            return destination_container;
        }

        template<typename key_type,typename value_type, typename source_container_type>
        static void to_map(const source_container_type& source_container,
                           const shared_ptr<map<key_type,value_type>>& destination_container){
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container, destination_container->end()));
        }

        template<typename value_type, typename source_container_type>
        static shared_ptr<set<value_type>> to_set(source_container_type& source_container){
            auto destination_container = make_shared<set<value_type>>();
            destination_container->reserve(source_container->size());
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container,destination_container->end()));
            return destination_container;
        }

        template<typename value_type, typename source_container_type>
        static void to_set(const source_container_type& source_container,
                           const shared_ptr<set<value_type>>& destination_container){
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container,destination_container->end()));
        }

        template<typename key_type,typename value_type, typename source_container_type>
        static shared_ptr<unordered_map<key_type,value_type>> to_unordered_map(const source_container_type& source_container){
            auto destination_container = make_shared<unordered_map<key_type,value_type>>();
            destination_container->reserve(source_container->size());
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container,destination_container->end()));
            return destination_container;
        }

        template<typename key_type,typename value_type, typename source_container_type>
        static void to_unordered_map(const source_container_type& source_container,
                                     const shared_ptr<unordered_map<key_type,value_type>>& destination_container){
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container,destination_container->end()));
        }

        template<typename value_type, typename source_container_type>
        static shared_ptr<unordered_set<value_type>> to_unordered_set(const source_container_type& source_container){
            auto destination_container = make_shared<unordered_set<value_type>>();
            destination_container->reserve(source_container->size());
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container,destination_container->end()));
            return destination_container;
        }

        template<typename value_type, typename source_container_type>
        static void to_unordered_set(const source_container_type& source_container,
                                     const shared_ptr<unordered_set<value_type>>& destination_container){
            copy(source_container->begin(),source_container->end(),
                 inserter(*destination_container,destination_container->end()));
        }

        template<typename value_type, typename source_container_type>
        static shared_ptr<vector<value_type>> to_vector(const source_container_type& source_container){
            auto destination_container = make_shared<vector<value_type>>();
            destination_container->reserve(source_container->size());
            copy(source_container->begin(),source_container->end(),
                 back_inserter(*destination_container));
            return destination_container;
        }

        template<typename value_type, typename source_container_type>
        static void to_vector(const source_container_type& source_container,
                              const shared_ptr<vector<value_type>>& destination_container){
            copy(source_container->begin(),source_container->end(),
                 back_inserter(*destination_container));
        }
    };
}