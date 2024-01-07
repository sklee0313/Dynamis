#include "PreProcessing.h"

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