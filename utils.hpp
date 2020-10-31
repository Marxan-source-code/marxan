#pragma once
// Helper functions for marxan execution.

#include <algorithm>
#include <cctype>
#include <locale>
#include <random>

namespace marxan {
namespace utils {

// return random number between 0 and parameter n-1
// not recommended to call repeatedly as the initialization is relatively expensive. Use for one-off calls.
// if needed for repeated use, consider initializing the distribution locally.
inline
int returnRandom (int num, std::mt19937& mt)
{
    std::uniform_int_distribution<int> int_range(0, num-1);
    return int_range(mt);
}

// similar to above, not recommended for repeated use in loops.
inline
double returnRandomFloat(std::mt19937& mt) {
    std::uniform_real_distribution<double> float_range(0.0, 1.0);
    return float_range(mt);
}

inline
string getFileNameSuffix(int suffixMode) {
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
inline
void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline
void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
inline
void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// Adds '/' or '\' to end of dir string if not existent
inline
string cleanDirectoryString(string dirName) {
    if (dirName.back() != '\\' && dirName.back() != '/') {
        return dirName + "/";
    }

    return dirName;
}

} // namespace utils
} // namespace marxan