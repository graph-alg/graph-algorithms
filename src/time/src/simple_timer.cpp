
#include "time/simple_timer.h"

namespace scnu
{
    /**
     * @brief construct a simple_timer
     */
    simple_timer::simple_timer()
    {
        start_time = high_resolution_clock::now();
    }

    /**
     * @fn get_elapse_millisecond
     * @brief get a millisecond result
     * @return
     */
    double simple_timer::get_elapse_millisecond()
    {
        end_time = high_resolution_clock::now();
        return duration<double,ratio<1,1000>>(end_time - start_time).count();
    }

    /**
     * @fn get_elapse_second
     * @brief get a second result
     * @return
     */
    double simple_timer::get_elapse_second()
    {
        end_time = high_resolution_clock::now();
        return duration<double,ratio<1,1>>(end_time - start_time).count();
    }
}