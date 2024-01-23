/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : wing_compare.h
* @brief      : common header files for k-core
* @version    : 1.0
* @date       : 2022/5/10
******************************************************************************************************************/

#include "wing/wing_utility.h"

namespace scnu{
    class wing_compare {
    public:
        static bool same_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container1,
                                         const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container2);

        static bool sub_associative_map(const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container1,
                                        const shared_ptr<unordered_map<shared_ptr<abstract_bipartite_edge>, uint32_t>>& container2);
    };
}



