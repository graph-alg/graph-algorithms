/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : string_algorithm.h
* @brief      : Some useful functions for string
* @version    : 1.0
* @date       : 2019/09/02
******************************************************************************************************************/

#pragma once
#include "string/string_utility.h"

namespace scnu
{
    /**
     * @class string_algorithm
     * @brief some useful functions for string
     */
    class string_algorithm{
    public:
        template< typename... Args >
        string format(const char* format, Args... args) {
            uint32_t length = std::snprintf(nullptr, 0, format, args...);
            if (length <= 0) {
                return "";
            }

            auto buf = std::make_unique<char>(length + 1);
            std::snprintf(buf.get(), length + 1, format, args...);

            return string(buf.get());
        }

        static shared_ptr<vector<string>> find_split(const string& original_string,
                                                                    const string& delimiter);

        static string left_trim(const string &original_string, char c);

        static shared_ptr<vector<string>> regex_split(const string& original_string,
                                                                     const string& delimiter);

        static string replace_all(const string& original_string, const string& old_substring,
                                  const string& new_substring);

        static string right_trim(const string &original_string, char c);

        static shared_ptr<vector<string>> strtok_split(const string& original_string,
                                                                      const string& delimiter);

        static shared_ptr<vector<string>> stream_split(const string& original_string,
                                                                      const string& delimiter);

        static string trim(const string &original_string, char c);

    };
}


