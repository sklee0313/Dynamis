#include <stdexcept>
#include <iomanip>
#include "PreProcessing.h"
#include "LinearElasticity.h"
using namespace Dynamis::PreProcessing;
namespace Dynamis::ConstitutiveLaw
{
    void LinearElasticity::computeC()
    {
        Eigen::Matrix<double, 3, 3> C11;
        Eigen::Matrix<double, 3, 3> C22;

        // Compute C11
        C11 << 1.0, nu / (1.0 - nu), nu / (1.0 - nu),
            nu / (1.0 - nu), 1.0, nu / (1 - nu),
            nu / (1.0 - nu), nu / (1.0 - nu), 1.0;

        // Compute C22
        C22 = Eigen::Matrix<double, 3, 3>::Identity();
        C22 *= (1 - 2 * nu) / (2 - 2 * nu);

        // Assemble matrix C
        C.topLeftCorner(3, 3) = C11;
        C.bottomRightCorner(3, 3) = C22;
        C *= E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    }

    LinearElasticity::LinearElasticity(std::unique_ptr<std::ifstream> &file)
    {
        // Show we use linear elasticity
        std::cout << GREEN << std::left << std::setw(20)
                  << "Constitutive law: "
                  << "Linear Elasticity" << RESET << std::endl;

        // Read nu, E, and density
        bool nuRead = false, ERead = false, densityRead = false;

        std::string line;
        (*file).seekg(0, std::ios::beg); // Reset to the beginning of the file
        while (std::getline(*file, line) && !(nuRead && ERead && densityRead))
        {
            if (!ERead && Dynamis::PreProcessing::tryReadValue<double>(file, line, "YoungsModulus", E))
            {
                ERead = true;
                std::cout << GREEN << std::left << std::setw(20)
                          << "Young's modulus: " << E / 1e9 << " GPa" << RESET << std::endl;
            }
            else if (!nuRead && Dynamis::PreProcessing::tryReadValue<double>(file, line, "PoissonsRatio", nu))
            {
                nuRead = true;
                std::cout << GREEN << std::left << std::setw(20)
                          << "Poisson's ratio: " << nu << RESET << std::endl;
            }
            else if (!densityRead && Dynamis::PreProcessing::tryReadValue<double>(file, line, "Density", density))
            {
                densityRead = true;
                std::cout << GREEN << std::left << std::setw(20)
                          << "Density: " << density << " kg/m^3" << RESET << std::endl;
            }
        }

        if (!nuRead || !ERead || !densityRead) // if
        {
            std::string errorMessage = RED + "Failed to read necessary material properties: ";

            if (!nuRead)
            {
                errorMessage += "Poisson's ratio ";
            }
            if (!ERead)
            {
                errorMessage += "Young's modulus ";
            }
            if (!densityRead)
            {
                errorMessage += "Density ";
            }

            throw std::runtime_error(errorMessage);
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