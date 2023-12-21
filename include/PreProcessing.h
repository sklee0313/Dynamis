#ifndef PreProcessing_H
#define PreProcessing_H

#include <fstream>
#include <iostream>
#include <cstdlib> // Include for std::exit and EXIT_FAILURE
#include <algorithm>
#include <cctype>
#include <locale>

namespace Dynamis::PreProcessing
{
    // ANSI Escape Code for Colors
    const std::string RED = "\033[31m";   // Red
    const std::string GREEN = "\033[32m"; // Green
    const std::string RESET = "\033[0m";  // Reset to default color

    // Trim from start (in place)
    static inline void ltrim(std::string &s)
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
                                        { return !std::isspace(ch); }));
    }

    // Trim from end (in place)
    static inline void rtrim(std::string &s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch)
                             { return !std::isspace(ch); })
                    .base(),
                s.end());
    }

    // Trim from both ends (in place)
    static inline void trim(std::string &s)
    {
        ltrim(s);
        rtrim(s);
    }

    bool tryReadValue(std::ifstream &file, std::string &line, const std::string &identifier, double &value);

    std::ifstream InputChecker(int argc, char *argv[]);
}

#endif // PreProcessing_H