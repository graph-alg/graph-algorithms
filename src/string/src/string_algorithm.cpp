/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : extend_string.cpp
* @brief      : A subclass of string for implementing extra functions
* @version    : 1.1
* @date       : 2019/10/09
******************************************************************************************************************/

#include "string/string_algorithm.h"

namespace scnu {

    /**
     * @brief split the original_string based on the STL find method
     * @param original_string
     * @param delimiter
     * @return
     */
    shared_ptr<vector<string>> string_algorithm::find_split(const string &original_string,
                                                            const string &delimiter) {
        auto string_vector = make_shared<vector<string>>();
        uint32_t begin = 0;
        while (begin < original_string.length()) {
            uint32_t end = original_string.find(delimiter, begin);
            if (end < original_string.length()) {
                string substring = move(original_string.substr(begin, (end - begin)));
                string_vector->emplace_back(substring);
                begin = end + 1;
            } else {
                string substring = move(original_string.substr(begin, (end - begin)));
                string_vector->emplace_back(substring);
                break;
            }
        }
        return string_vector;
    }

    /**
     * @brief trim the give character from the string at left side
     * @param original_string
     * @param c
     * @return
     */
    string string_algorithm::left_trim(const string &original_string, char c) {
        auto index = original_string.find_first_of(c);
        string substr{original_string};
        while (index != string::npos) {
            substr = original_string.substr(index + 1);
            index = substr.find_first_of(c);
        }
        return substr;
    }

    /**
     * @brief split the original string based on the regex
     * @param original_string
     * @param delimiter
     * @return
     */
    shared_ptr<vector<string>> string_algorithm::regex_split(const string &original_string,
                                                             const string &delimiter) {
        regex substring_regex{delimiter};
        return make_shared<vector<string>>(
                sregex_token_iterator(original_string.begin(), original_string.end(), substring_regex, -1),
                sregex_token_iterator()
        );
    }

    /**
     * @brief replace all old strings with the new string
     * @note it does not revise the original string
     * @param original_string
     * @param old_substring
     * @param new_substring
     * @return
     */
    string string_algorithm::replace_all(const string &original_string, const string &old_substring,
                                         const string &new_substring) {
        regex string_regex(old_substring);
        return regex_replace(original_string, string_regex, new_substring);
    }

    /**
     * @brief split the original string based on a C-style function strtok
     * @param original_string
     * @param delimiter
     * @return
     */
    shared_ptr<vector<string>> string_algorithm::strtok_split(const string &original_string,
                                                              const string &delimiter) {
        auto string_vector = make_shared<vector<string>>();
        auto string_array = make_shared<char>(original_string.length() + 1);
        strcpy(string_array.get(), original_string.c_str());
        const char *delimiter_pointer = delimiter.c_str();
        const char *p = strtok(string_array.get(), delimiter_pointer);
        while (p) {
            string substring = p;
            string_vector->emplace_back(substring);
            p = strtok(nullptr, delimiter_pointer);
        }
        return string_vector;
    }

    /**
     * @brief split the original string based on string stream
     * @param original_string
     * @param delimiter
     * @return
     */
    shared_ptr<vector<string>> string_algorithm::stream_split(const string &original_string,
                                                              const string &delimiter) {
        auto string_vector = make_shared<vector<string>>();
        stringstream string_stream(original_string);
        string substring;
        while (getline(string_stream, substring, *delimiter.c_str())) {
            string_vector->emplace_back(substring);
        }
        return string_vector;
    }


    /**
     * @brief trim the give character from the string at right side
     * @param original_string
     * @param c
     * @return
     */
    string string_algorithm::right_trim(const string &original_string, char c) {
        auto index = original_string.find_last_of(c);
        string substr{original_string};
        while (index != string::npos) {
            substr = original_string.substr(0, index - 1);
            index = substr.find_last_of(c);
        }
        return substr;
    }

    /**
     * @brief trim the give character from the string at both sides
     * @param original_string
     * @param c
     * @return
     */
    string string_algorithm::trim(const string &original_string, char c) {
        auto substr{right_trim(original_string, c)};
        return left_trim(substr, c);
    }
}
