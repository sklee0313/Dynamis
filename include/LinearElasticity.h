#ifndef LINEAR_ELASTICITY_H
#define LINEAR_ELASTICITY_H

#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <iostream>

namespace Dynamis::ConstitutiveLaw
{

    class LinearElasticity
    {
    private:
        double nu, E, density;
        Eigen::Matrix<double, 6, 6> C;

        void computeC();

    public:
        LinearElasticity(std::ifstream &file);
        const Eigen::Matrix<double, 6, 6> &getC() const;
        const double &getDensity() const;
    };

    // More laws to come..

}

#endif // LINEAR_ELASTICITY_H
