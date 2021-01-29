#pragma once
// Helper functions for marxan execution.
// Unit functions with no depedency on marxan specific structures or logic.

#include <algorithm>
#include <cctype>
#include <locale>
#include <random>
#include <string>
#include <sstream>

namespace marxan {
    namespace utils {

        inline std::string getFileNameSuffix(int suffixMode) {
            if (suffixMode == 3) {
                return ".csv";
            }
            else if (suffixMode == 2) {
                return ".txt";
            }
            else {
                return ".dat";
            }
        }

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


    } // namespace utils
} // namespace marxan