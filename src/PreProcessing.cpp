#include "PreProcessing.h"

bool Dynamis::PreProcessing::tryReadMatrix(std::unique_ptr<std::ifstream> &file, std::string &line, const std::string &identifier, Eigen::MatrixXd &matrix, const size_t &rows, const size_t &cols)
{
    std::string openingTag = "<" + identifier + ">";
    std::string closingTag = "</" + identifier + ">";

    trim(line);
    std::string token;
    if (openingTag == line)
    {
        matrix.resize(rows, cols);
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                if (!(*file >> matrix(i, j)))
                {
                    std::string errorMessage = "Error with reading matrix(" + std::to_string(i) + ", " + std::to_string(j) + ")";
                    throw std::runtime_error(errorMessage);
                };
            }
        }
        *file >> token;
        if (token == closingTag)
        {
            return true; // Successfully read the value
        }
    }
    return false;
}

std::unique_ptr<std::ifstream> Dynamis::PreProcessing::InputChecker(int argc, char *argv[])
{

    if (argc != 2)
    {
        throw std::invalid_argument("Usage: " + std::string(argv[0]) + " <path_to_mesh_file>");
    }
    auto file = std::make_unique<std::ifstream>(argv[1]);

    if (!file->is_open())
    {
        throw std::runtime_error("Input file " + std::string(argv[1]) + " cannot be found");
    }
    return file; // Return the ifstream object if no exception was thrown
}

bool Dynamis::PreProcessing::readMatrix(std::unique_ptr<std::ifstream> &file, const std::string &key, Eigen::MatrixXd &matrix, const size_t &numRows, const size_t &numCols)
{
    std::string line;
    bool matrixRead = false;
    (*file).seekg(0, std::ios::beg); // Reset to the beginning of the file
    while (std::getline(*file, line) && !matrixRead)
    {
        if (Dynamis::PreProcessing::tryReadMatrix(file, line, key, matrix, numRows, numCols))
        {
            matrixRead = true;
            std::cout << "Successful to read " + key << std::endl;
            break;
        }
    }
    return matrixRead;
}