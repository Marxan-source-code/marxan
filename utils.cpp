// Helper functions for marxan execution.
// Unit functions with no depedency on marxan specific structures or logic.

#include <algorithm>
#include <cctype>
#include <locale>
#include <string>
#include <vector>
#include <sstream>
#include "utils.hpp"

namespace marxan {
    namespace utils {

        std::string getFileNameSuffix(int suffixMode) {
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

        std::vector<std::string> get_tokens(const std::string& str)
        {
            static const std::string delimeters(",;\t");
            std::vector<std::string> tokens;
            std::string word;
            for (char ch : str)
            {
                if (delimeters.find_first_of(ch) == std::string::npos)
                    word.push_back(ch);
                else
                {
                    if (!word.empty())
                    {
                        tokens.push_back(word);
                        word.clear();
                    }
                }
            }
            if (!word.empty())
                tokens.push_back(word);
            return tokens;
        }

        std::stringstream stream_line(const std::string& str)
        {
            static const std::string delimeters(" ,;:^*\"/\t\'\\\n");
            std::stringstream ss;
            for (char ch : str)
            {
                if (delimeters.find_first_of(ch) == std::string::npos)
                    ss << ch;
                else
                {
                    ss << ' ';
                }
            }
            return ss;
        }

        bool is_like_numerical_data(const std::string& str)
        {
            //check if there are chars other then delimeters and ones used to represent numbers
            size_t letter_pos = str.find_first_not_of(" ,;:^*\"/\t\'\\\n\r.+-Ee0123456789");
            if (letter_pos != std::string::npos)
                return false;
            return true;
        }

        char guess_delimeter(const std::string& str)
        {
            char delim = ',';
            int n_commas = std::count(str.begin(), str.end(), ',');
            int n_tabs = std::count(str.begin(), str.end(), '\t');
            int n_semicolons = std::count(str.begin(), str.end(), ';');
            if(n_tabs > n_commas)
                delim = '\t';
            if(n_semicolons > n_commas && n_semicolons > n_tabs)
                delim = ';';
            return delim;
        }


    } // namespace utils
} // namespace marxan
