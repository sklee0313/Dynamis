#include "PreProcessing.h"

bool Dynamis::PreProcessing::tryReadMatrix(std::ifstream &file, std::string &line, const std::string &identifier, Eigen::MatrixXd &matrix, const size_t &rows, const size_t &cols)
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
                if (!(file >> matrix(i, j)))
                {
                    std::string errorMessage = "Error with reading matrix(" + std::to_string(i) + ", " + std::to_string(j) + ")";
                    throw std::runtime_error(errorMessage);
                };
            }
        }
        file >> token;
        if (token == closingTag)
        {
            return true; // Successfully read the value
        }
    }
    return false;
}

std::ifstream Dynamis::PreProcessing::InputChecker(int argc, char *argv[])
{
    std::ifstream infile; // Declare infile outside the try block to return it

    try
    {
        if (argc != 2)
        {
            throw std::invalid_argument("Usage: " + std::string(argv[0]) + " <path_to_mesh_file>");
        }

        infile.open(argv[1]);

        if (!infile)
        {
            throw std::runtime_error("Input file " + std::string(argv[1]) + " cannot be found");
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE); // Exit the program with a failure status
    }

    return infile; // Return the ifstream object if no exception was thrown
}