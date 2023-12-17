#include "LinearElasticity.hpp"
#include "Utility.hpp"
#include <stdexcept>

namespace EFEM::ConstitutiveLaw
{

    bool LinearElasticity::tryReadValue(std::ifstream &file, std::string &line, const std::string &identifier, double &value)
    {
        trim(line);
        if (identifier == line)
        {
            return bool(file >> value);
        }
        return false;
    }

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

    LinearElasticity::LinearElasticity(std::ifstream &file)
    {
        bool nuRead = false, ERead = false, densityRead = false;

        std::string line;
        while (std::getline(file, line) && !(nuRead && ERead && densityRead))
        {
            if (!ERead && tryReadValue(file, line, "YoungsModulus", E))
            {
                ERead = true;
                std::cout << ERead << std::endl;
            }
            else if (!nuRead && tryReadValue(file, line, "PoissonsRatio", nu))
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