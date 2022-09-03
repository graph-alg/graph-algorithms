/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : random_generator.h
* @brief      : A class for producing random integers
* @version    : 1.0
* @date       : 2019/09/02
******************************************************************************************************************/

#pragma once
#include "random/random_utility.h"

namespace scnu
{
    /**
     * @class random_generator
     * @brief a class for random number generator
     */
    class random_generator {
    public:
        random_generator();

        explicit random_generator(const shared_ptr<random_device>& other_random_device);

        template<typename engine_type>
        explicit random_generator(const shared_ptr<random_device>& other_random_device, const shared_ptr<engine_type>& other_engine):
                default_random_device(other_random_device),default_engine(other_engine)
        {

        }

        static shared_ptr<default_random_engine> get_default_engine(const shared_ptr<random_device>& rd);

        static shared_ptr<mt19937> get_mt19937_engine(const shared_ptr<random_device>& rd);

        int get_random_integer(const shared_ptr<uniform_int_distribution<int>>& distribution);

        static shared_ptr<uniform_int_distribution<int>> get_uniform_distribution(int min_number, int max_number);

        shared_ptr<normal_distribution<double>> get_normal_distribution(double mu, double sigma);

    private:
        shared_ptr<random_device> default_random_device;
        shared_ptr<default_random_engine> default_engine;
    };
}



