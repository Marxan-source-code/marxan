#pragma once
// Helper functions for marxan execution.
// Unit functions with no depedency on marxan specific structures or logic.

#include <algorithm>
#include <cctype>
#include <locale>
#include <random>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

namespace marxan {
    namespace utils {

        std::string getFileNameSuffix(int suffixMode);

        // String trimming functions credited to https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
        // trim from start (in place)
        inline void ltrim(std::string& s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
                return !std::isspace(ch);
                }));
        }

        // trim from end (in place)
        inline void rtrim(std::string& s) {
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
                }).base(), s.end());
        }

        // trim from both ends (in place)
        inline void trim(std::string& s) {
            ltrim(s);
            rtrim(s);
        }

        inline void to_lower(std::string& str)
        {
            transform(str.begin(), str.end(), str.begin(), ::tolower);
        }

        // converts number to a string padded with leading zeros
        // does nothing if stringLength is less than the digits in number.
        inline std::string intToPaddedString(int number, int stringLength)
        {
            std::ostringstream ss;
            ss << std::setw(stringLength) << std::setfill('0') << number;
            return ss.str();
        }

        // Adds '/' to end of dir string if not existent
        inline std::string cleanDirectoryString(std::string dirName) {
            if (dirName.back() != '\\' && dirName.back() != '/') {
                return dirName + "/";
            }

            return dirName;
        }

        inline double probZUT(double z)
            /*
            Using erf as it can be transformed into the standard normal function:
            http://www.cplusplus.com/reference/cmath/erf/
            See discussion https://stackoverflow.com/questions/27214780/how-to-implement-the-standard-normal-cumulative-distribution-function-in-c-or-o
            */
        {
            return 1 - (1 / 2) * (1 + erf(z / sqrt(2)));
        }

        inline double probZLT(double z)
        {
            return (1 / 2) * (1 + erf(z / sqrt(2)));
        }

        // given the parsed value of an option and its given value storage, parse it in.
        // For string types, the overload below is given.
        inline
            void readInputOptionValue(std::stringstream& parsed, std::string& value) {
            value = parsed.str();
        }

        template<class T> inline
            void readInputOptionValue(std::stringstream& parsed, T& value) {
            parsed >> value;
        }

        class formatted_string_stream
        {
            public:
            formatted_string_stream(const std::string& s, char delim)
            {
                s_ = s;
                pos_ = 0;
                delim_ = delim;
                fail_ = false;
                
            }

            bool fail()
            {
                return fail_;
            }

            formatted_string_stream& operator>>(std::string& s) 
            {
                if(pos_ >= s_.size())
                {
                    fail_ = true;
                    return *this;
                }

                size_t end_token = this->token_end();
                s.clear();
                for(size_t i = pos_; i < end_token; i++ )
                    if(s_[i] != '"')
                        s.push_back(s_[i]);
                utils::trim(s);
                size_t next_pos = end_token + 1; //skip delimeter
                pos_ = std::min(next_pos, s_.size());
                return *this;                  
            }

            
            formatted_string_stream& operator>>(double& val) 
            {
                size_t end_token = token_end();
                try
                {
                    val = stod(s_.substr(pos_, end_token));
                }
                catch(const std::invalid_argument& e)
                {
                    fail_ = true;    
                }
                size_t next_pos = end_token + 1; //skip delimeter
                pos_ = std::min(next_pos, s_.size());
                return *this;
            }

            formatted_string_stream& operator>>(int& n) 
            {
                size_t end_token = token_end();
                try
                {
                    n = stoi(s_.substr(pos_, end_token));
                }
                catch(const std::invalid_argument& e)
                {
                        fail_ = true;    
                }
                size_t next_pos = end_token + 1; //skip delimeter
                pos_ = std::min(next_pos, s_.size());
                return *this;
            }

            private:

            size_t token_end()
            {
                size_t next_pos = s_.find_first_of(delim_, pos_);
                if(next_pos == std::string::npos)
                    next_pos = s_.size();
                return next_pos;
            }

            std::string s_;
            size_t pos_;
            char delim_;
            bool fail_;
        };

        //Split string on tokens using delimeters
        std::vector<std::string> get_tokens(const std::string& str, char delim);

        //check if string likely repesents only some numbers
        bool is_like_numerical_data(const std::string& str);

        //check if string likely repesents only some numbers
        char guess_delimeter(const std::string& str);


    } // namespace utils
} // namespace marxan