#include "ConstitutiveLaw.hpp"

namespace EFEM::ConstitutiveLaw
{
    Eigen::Matrix<double, 6, 6> LinearElasticity(std::ifstream &infile)
    {
        double nu, E;
        infile >> nu >> E;

        Eigen::Matrix<double, 6, 6> C;
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

        return C;
    }
}
