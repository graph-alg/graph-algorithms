//
// Created by wbai on 22-5-18.
//

#include "wing/wing_compare.h"

namespace scnu{
    bool wing_compare::same_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container1,
                                            const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container2)
    {
        if(container1->size()!=container2->size()){
            printf("error size!");
            return false;
        }
        return sub_associative_map(container1, container2)
               && sub_associative_map(container2, container1);
    }

    bool wing_compare::sub_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container1,
                                           const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container2)
    {
        return std::all_of(begin(*container1), end(*container1),
                           [&](const auto &iter) {
                               if (!container2->count(iter.first)|| container2->at(iter.first) != iter.second) {
                                   printf("%u,%u,%u,%u\n",iter.first->get_left_vertex_id(),iter.first->get_right_vertex_id() ,iter.second,container2->at(iter.first));
                                   return false;
                               }
                               return true;
                           });
    }
}
