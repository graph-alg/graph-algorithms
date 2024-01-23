/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : function_wrapper.h
* @brief      : A class for wrapper function
* @version    : 1.0
* @date       : 2020/9/2
******************************************************************************************************************/

#pragma once
#include "thread/thread_utility.h"

namespace scnu
{
    /**
     * @class function_wrapper
     * @brief A class for wrapper function
     */
    class function_wrapper
    {
    public:
        template<typename function_type>
        explicit function_wrapper(function_type&& function)
                :impl(new impl_type<function_type>(std::forward<function_type>(function)))
        {

        }

        void operator()()
        {
            impl->call();
        }

        function_wrapper()=default;

        function_wrapper(function_wrapper&& function_wrapper)noexcept
                  : impl(std::move(function_wrapper.impl))
        {

        }

        function_wrapper& operator=(function_wrapper&& function_wrapper) noexcept
        {
            impl = std::move(function_wrapper.impl);
            return *this;
        }

        function_wrapper(const function_wrapper&) = delete;
        function_wrapper(function_wrapper &) = delete;
        function_wrapper& operator=(const function_wrapper&) = delete;
    private:
        struct impl_base
        {
            virtual void call()=0;
            virtual ~impl_base() = default;
        };
        std::unique_ptr<impl_base> impl;

        template<typename  function_type>
        struct impl_type: impl_base
        {
            function_type function;
            explicit impl_type(function_type&& other_function): function(std::move(other_function)){

            }
            void call() override{
                function();
            };
        };
    };
}




