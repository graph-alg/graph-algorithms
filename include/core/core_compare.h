/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : core_compare.h
* @brief      : common header files for k-core
* @version    : 1.0
* @date       : 2022/5/10
******************************************************************************************************************/

#include "core/core_utility.h"

namespace scnu{
    class core_compare {
    public:
        static bool same_associative_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& container1,
                                         const shared_ptr<unordered_map<uint32_t, uint32_t>>& container2);

        static bool sub_associative_map(const shared_ptr<unordered_map<uint32_t, uint32_t>>& container1, const shared_ptr<unordered_map<uint32_t, uint32_t>>& container2);
    };
}


