#ifndef LINEAR_ELASTICITY_HPP
#define LINEAR_ELASTICITY_HPP

#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <iostream>

namespace EFEM::ConstitutiveLaw
{

    class LinearElasticity
    {
    private:
        double nu, E, density;
        Eigen::Matrix<double, 6, 6> C;

        bool tryReadValue(std::ifstream &file, std::string &line, const std::string &identifier, double &value);
        void computeC();

    public:
        LinearElasticity(std::ifstream &file);
        const Eigen::Matrix<double, 6, 6> &getC() const;
        const double &getDensity() const;
    };

}

#endif // LINEAR_ELASTICITY_HPP
