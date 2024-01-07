#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "Element.h"
namespace Dynamis
{
    using namespace Eigen;

    class OFE8nodeLinear
    {

    private:
        double delh = 5e-5;
        double beta = 0.01;
        double density;
        Matrix3d J, invJ;

        MatrixXd dHdR, dHdx, dHdy, dHdz, weight;
        Matrix<double, 6, 96> B;
        Matrix<double, 96, 96> K;
        Matrix<double, 32, 32> M;
        Matrix<double, 1, 32> H;
        Matrix<double, 1, 8> rho_p, rho_m;
        MatrixXd rho_r, rho_s, rho_t, rho, h, h_r, h_s, h_t;

        double r, s, t, x, y, z, detJ, rhoK_x, rhoK_y, rhoK_z, xK, yK, zK;
        int row, col;

    public:
        OFE8nodeLinear(const MatrixXd &NQ);
        int nodaldofs = 12;
        void ElementStiffness(const MatrixXd &C, const MatrixXd &nodes, const Element &element, std::vector<Triplet<double>> &triplets, const MatrixXd &NQ, const std::vector<double> &radius);
        void ElementMass(const double &density, const MatrixXd &nodes, const Element &element, std::vector<Triplet<double>> &triplets, const MatrixXd &NQ, const std::vector<double> &radius);
        void CalRho(const double &r, const double &s, const double &t, Matrix<double, 1, 8> &rho);
        void Jacobian(const MatrixXd &nodes, const Element &element, const int &i);
        void OFEintp(Matrix<double, 1, 32> &H, const double &x, const double &y, const double &z, const int &k, const MatrixXd &nodes, const Element &element, const std::vector<double> &radius);
        std::vector<double> radius(const MatrixXd &nodes, const std::vector<Element> &elements, const size_t &ne, const size_t &nn);
    };

}