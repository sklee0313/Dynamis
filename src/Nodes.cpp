#include <stdexcept>

#include "Nodes.h"
#include "PreProcessing.h"
namespace Dynamis::core
{
    Nodes::Nodes(std::ifstream &file)
    {
        bool nuRead = false, ERead = false, densityRead = false;

        std::string line;
        file.seekg(0, std::ios::beg); // Reset to the beginning of the file
        while (std::getline(file, line) && !(nuRead && ERead && densityRead))
        {
            if (!ERead && Dynamis::PreProcessing::tryReadValue(file, line, "Number of Nodes", E))
            {
                ERead = true;
                std::cout << ERead << std::endl;
            }
            else if (!nuRead && Dynamis::PreProcessing::tryReadValue(file, line, "Nodes", nu))
            {
                nuRead = true;
                std::cout << nu << std::endl;
            }
            else if (!densityRead && tryReadValue(file, line, "Density", density))
            {
                densityRead = true;
                std::cout << density << std::endl;
            }
        }

        if (!nuRead || !ERead || !densityRead)
        {
            std::cout << nuRead << ERead << densityRead << std::endl;
            std::cout << nu << " " << E << " " << density << std::endl;
            throw std::runtime_error("Failed to read necessary material properties from the file");
        }

        computeC();
    }

    const Eigen::Matrix<double, 6, 6> &LinearElasticity::getC() const
    {
        return C;
    }
    const double &LinearElasticity::getDensity() const
    {
        return density;
    }
}