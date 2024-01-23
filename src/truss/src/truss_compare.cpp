
#include "truss/truss_compare.h"

namespace scnu{
    bool truss_compare::same_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container1,
                                     const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container2)
    {
        if(container1->size()!=container2->size()){
            printf("error size!");
            return false;
        }
        return sub_associative_map(container1, container2)
               && sub_associative_map(container2, container1);
    }

    bool truss_compare::sub_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container1,
                                    const shared_ptr<unordered_map<shared_ptr<abstract_edge>, uint32_t>>& container2)
    {

        for(const auto &[e, e_truss_number]:*container1){
            if(!container2->count(e)){
                printf("Unequal Key: %u,%u\n",e->get_source_vertex_id(), e->get_destination_vertex_id());
                return false;
            }
            if(container2->at(e) != e_truss_number){
                printf("%u,%u:%u,%u\n",e->get_source_vertex_id(), e->get_destination_vertex_id(),e_truss_number,container2->at(e));
                return false;
            }
        }
        return true;
    }
}