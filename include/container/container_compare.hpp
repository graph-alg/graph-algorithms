/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : HashMapCompare.hpp
* @brief      : A template class for comparing elements in given containers
* @version    : 1.1
* @date       : 2020/10/12
******************************************************************************************************************/

#pragma once
#include "container/container_utility.h"

namespace scnu
{
    /**
     * @details common functions for comparing elements in containers
     * @remarks elements in two containers should have the sub type and operator =
     */
    class container_compare
    {
    public:
        /**
         * @brief compare the elements of two map are equal or not
         * @note parameters should be smart point type
         * @tparam T
         * @param container1
         * @param container2
         * @return
         */
        template<typename key_type, typename value_type>
        static bool same_associative_map(const shared_ptr<unordered_map<key_type,value_type>>& container1,
                                         const shared_ptr<unordered_map<key_type,value_type>>& container2)
        {
            if(container1->size()!=container2->size()){
                printf("error size!");
                return false;
            }
            return sub_associative_map(container1, container2)
                   && sub_associative_map(container2, container1);
        }

        /**
         * @brief compare the elements of two associated container are equal or not
         * @note parameters should be smart point type
         * @tparam T
         * @param container1
         * @param container2
         * @return
         */
        template<typename T>
        static bool same_associative_set(T container1, T container2)
        {
            if(container1->size()!=container2->size()){
                printf("error size!");
                return false;
            }
            return sub_associative_set(container1, container2)
                   && sub_associative_set(container2, container1);
        }


        /**
         * @brief compare the elements of two sequence container are equal or not
         * @note parameters should be smart point type
         * @tparam T
         * @param container1
         * @param container2
         * @return
         */
        template<typename T>
        static bool same_sequence_container(T container1, T container2)
        {
            return sub_sequence_container(container1, container2)
                   && sub_sequence_container(container2, container1);
        }

        /**
         * @brief check the element of one associative map is a subset of the other
         * @note parameters should be smart pointer type
         * @tparam T
         * @param container1
         * @param container2
         * @return
         */
        template<typename T>
        static bool sub_associative_map(const T& container1, const T& container2)
        {
            return std::all_of(begin(*container1), end(*container1),
                               [&](const auto &iter) {
                                   if (!container2->count(iter.first)|| container2->at(iter.first) != iter.second) {
                                       printf("%u,%u,%u\n",iter.first,iter.second,container2->at(iter.first));
                                       return false;
                                   }
                                   return true;
                               });
        }

        /**
         * @brief check the element of one associative set is a subset of the other
         * @note parameters should be smart pointer type
         * @tparam T
         * @param container1
         * @param container2
         * @return
         */
        template<typename T>
        static bool sub_associative_set(const T& container1, const T& container2)
        {
            return std::all_of(begin(*container1), end(*container1),
                               [&](const auto &iter) {
                                   if (!container2->count(iter)) {
                                       printf("Unequal Key\n");
                                       return false;
                                   }
                                   return true;
                               });
        }

        /**
         * @brief check the element of one sequence container is a subset of the other
         * @note parameters should be smart pointer type
         * @tparam T
         * @param container1
         * @param container2
         * @return
         */
        template<typename T>
        static bool sub_sequence_container(T container1, T container2)
        {
            return all_of(begin(*container1), end(*container1),
                               [&](const auto &element) {
                                   return container2->find(element);
                               });
        }
    };
}
