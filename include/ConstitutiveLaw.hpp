#ifndef CONSTITUTIVE_LAW_HPP
#define CONSTITUTIVE_LAW_HPP

#include <fstream>
#include <Eigen/Dense> // Include the Eigen library for Matrix operations

namespace EFEM::ConstitutiveLaw
{
    Eigen::Matrix<double, 6, 6> LinearElasticity(std::ifstream &infile);
}

#endif // CONSTITUTIVE_LAW_HPP
