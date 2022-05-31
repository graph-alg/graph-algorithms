//
// Created by wbai on 22-5-18.
//

#include "core/core_compare.h"

namespace scnu{
    bool core_compare::same_associative_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& container1,
                                     const shared_ptr<unordered_map<uint32_t, uint32_t>>& container2)
    {
        if(container1->size()!=container2->size()){
            printf("error size!");
            return false;
        }
        return sub_associative_map(container1, container2)
               && sub_associative_map(container2, container1);
    }

    bool core_compare::sub_associative_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& container1, const shared_ptr<unordered_map<uint32_t, uint32_t>>& container2)
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
}