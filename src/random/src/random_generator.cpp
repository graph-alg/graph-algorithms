/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : random_generator.cpp
* @brief      : An implementation of functions for random_generator
* @version    : 1.1
* @date       : 2019/10/08
******************************************************************************************************************/

#include "random/random_generator.h"

namespace scnu
{
    random_generator::random_generator() {
        default_random_device = make_shared<random_device>();
        default_engine = get_default_engine(default_random_device);
    }

    random_generator::random_generator(const shared_ptr<random_device>& other_random_device)
    :default_random_device(other_random_device)
    {

    }

    /**
     * @brief get a default default_engine
     * @return
     */
    shared_ptr<default_random_engine> random_generator::get_default_engine(const shared_ptr<random_device>& rd)
    {
        return make_shared<default_random_engine>((*rd)());
    }
    /**
     * @brief get a mt19937 default_engine for shuffle container
     * @return
     */
    shared_ptr<mt19937> random_generator::get_mt19937_engine(const shared_ptr<random_device>& rd)
    {
        return make_shared<mt19937>((*rd)());
    }

    /**
     * @brief a random number between with uniform distribution default_engine
     * @return
     */
    int random_generator::get_random_integer(const shared_ptr<uniform_int_distribution<int>>& distribution) {
        int number = (*distribution)(*default_engine);
        return number;
    }

    /**
     * @brief prepare for a uniform distribution integer generator
     * @param min_number
     * @param max_number
     */
    shared_ptr<uniform_int_distribution<int>> random_generator::get_uniform_distribution(int min_number, int max_number)
    {
        return make_shared<uniform_int_distribution<int>>(min_number, max_number);
    }

    /**
     * @brief prepare for a normal distribution double generator
     * @param mu
     * @param sigma
     * @return
     */
    shared_ptr<normal_distribution<double>> random_generator::get_normal_distribution(double mu, double sigma)
    {
        return make_shared<normal_distribution<double>>(mu,sigma);
    }
}