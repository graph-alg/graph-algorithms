/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : extend_string.h
* @brief      : A subclass of string for implementing extra functions
* @version    : 1.0
* @date       : 2022/03/18
******************************************************************************************************************/

#pragma once
#include "system_utility.h"


#define VMRSS_LINE 22 // the line number of memory
namespace scnu{
    class process_information {
    public:
        static uint32_t get_memory(pid_t p);
    };
}



