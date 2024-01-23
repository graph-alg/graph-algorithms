/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : extend_string.h
* @brief      : A subclass of string for implementing extra functions
* @version    : 1.0
* @date       : 2019/09/02
******************************************************************************************************************/

#pragma once
#include "string/string_utility.h"

namespace scnu
{
    /**
     * @class extend_string
     * @brief a subclass of string for implementing extra functions
     */
    class extend_string: public string {
        public:
            extend_string() = default;

            explicit extend_string(string  other_string);

            static string& get_delimiter();

            string& get_string();

            shared_ptr<vector<string>> split(const string& other_delimiter);

            friend istream& operator>>(istream& input_stream, extend_string& other_string);

            friend ostream& operator<<(ostream& output_stream, extend_string& other_string);
        private:
            static string delimiter;
            string data_string;
        };

        string extend_string::delimiter{" "};
}



