#ifndef PreProcessing_H
#define PreProcessing_H

#include <fstream>
#include <iostream>
#include <cstdlib> // Include for std::exit and EXIT_FAILURE
#include <algorithm>
#include <cctype>
#include <locale>
#include <memory>

#include "Eigen/Dense"

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

    bool tryReadMatrix(std::unique_ptr<std::ifstream> &file, std::string &line, const std::string &identifier, Eigen::MatrixXd &matrix, const size_t &rows, const size_t &cols);

    std::unique_ptr<std::ifstream> InputChecker(int argc, char *argv[]);

    template <typename T>
    bool tryReadValue(std::unique_ptr<std::ifstream> &file, std::string &line, const std::string &identifier, T &value)
    {
        std::string openingTag = "<" + identifier + ">";
        std::string closingTag = "</" + identifier + ">";

        trim(line);
        std::string token;
        if (openingTag == line)
        {
            if (*file >> value)
            {
                // Expect the closing tag
                *file >> token;
                if (token == closingTag)
                {
                    return true; // Successfully read the value
                }
            }
        }
        return false;
    }

    template <typename T>
    bool readValue(std::unique_ptr<std::ifstream> &file, const std::string &key, T &value)
    {
        std::string line;
        bool valueRead = false;
        (*file).seekg(0, std::ios::beg); // Reset to the beginning of the file
        while (std::getline(*file, line) && !valueRead)
        {
            if (Dynamis::PreProcessing::tryReadValue(file, line, key, value))
            {
                valueRead = true;
                std::cout << "Successful to read " << key << std::endl;
                break;
            }
        }
        return valueRead;
    }
    bool readMatrix(std::unique_ptr<std::ifstream> &file, const std::string &key, Eigen::MatrixXd &matrix, const size_t &numRows, const size_t &numCols);

}

#endif // PreProcessing_H