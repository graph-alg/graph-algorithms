/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : simple_timer.h
* @brief      : A simple simple_timer
* @version    : 1.0
* @date       : 2020/9/15
******************************************************************************************************************/

#pragma once
#include "time/time_utility.h"

namespace scnu
{
    /**
     * @class simple_timer
     * @brief a simple simple_timer for calculating used time
     */
    class simple_timer
    {
    public:
        simple_timer();

        double get_elapse_millisecond();

        double get_elapse_second();

    private:
        high_resolution_clock::time_point end_time;
        high_resolution_clock::time_point start_time;
    };
}




