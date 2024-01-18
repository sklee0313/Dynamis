#include "../include/OFE8nodeLinear.h"
#include "../include/Combination.h"
#include <iostream>
namespace Dynamis
{

    const Eigen::Matrix<double, 4, 2> OFE8nodeLinear::NQ_Mass = []
    {
        Eigen::Matrix<double, 4, 2> m;
        m << -0.861136311594052575224, 0.3478548451374538573731,
            -0.3399810435848562648027, 0.6521451548625461426269,
            0.3399810435848562648027, 0.6521451548625461426269,
            0.861136311594052575224, 0.3478548451374538573731;
        return m;
    }(); // no concern in extra copy: https://en.cppreference.com/w/cpp/language/copy_elision

    const Eigen::Matrix<double, 3, 2> OFE8nodeLinear::NQ_Stiffness = []
    {
        Eigen::Matrix<double, 3, 2> m;
        m << -0.7745966692414833770359, 0.5555555555555555555556,
            0, 0.8888888888888888888889,
            0.7745966692414833770359, 0.555555555555555555556;
        return m;
    }();

    void OFE8nodeLinear::initialization(const Eigen::MatrixXd &NQ)
    {
        B.setZero(); // strain-displacement matrix
        K.setZero(); // element stiffness matrix

        rho_r.resize(pow(NQ.rows(), 3.0), 8);  // derivative of PU w.r.t r-dir
        rho_s.resize(pow(NQ.rows(), 3.0), 8);  // derivative of PU w.r.t s-dir
        rho_t.resize(pow(NQ.rows(), 3.0), 8);  // derivative of PU w.r.t t-dir
        weight.resize(pow(NQ.rows(), 3.0), 1); // Gauss quadrature weights store

        rho.resize(pow(NQ.rows(), 3.0), 8); // PU functions store
        h.resize(pow(NQ.rows(), 3.0), 8);   // low-order FE function store
        h_r.resize(pow(NQ.rows(), 3.0), 8); // derivative of FE w.r.t r-dir
        h_s.resize(pow(NQ.rows(), 3.0), 8); // derivative of FE w.r.t s-dir
        h_t.resize(pow(NQ.rows(), 3.0), 8); // derivative of FE w.r.t. t-dir

        dHdx.resize(1, 32); // derivative of the displacement matrix w.r.t x-dir
        dHdy.resize(1, 32); // derivative of the displacement matrix w.r.t y-dir
        dHdz.resize(1, 32); // derivative of the displacement matrix w.r.t z-dir

        // Three FOR loops are used, each corresponds to a direction in three dimensions.
        // The iterator 'k' is the global iterator.
        // Numerical differentiation is used to evaluate derivatives of the PU functions.
        for (int i = 0, k = 0; i < NQ.rows(); i++)
        {
            for (int j = 0; j < NQ.rows(); j++)
            {
                for (int l = 0; l < NQ.rows(); l++, k++)
                {
                    weight(k) = NQ(i, 1) * NQ(j, 1) * NQ(l, 1);
                    // rho_r
                    s = NQ(j, 0);
                    t = NQ(l, 0);
                    r = NQ(i, 0) + delh;
                    CalRho(r, s, t, rho_p);
                    r = NQ(i, 0) - delh;
                    CalRho(r, s, t, rho_m);
                    rho_r(k, Eigen::all) = (rho_p - rho_m) / (2 * delh);

                    // rho_s
                    r = NQ(i, 0);
                    t = NQ(l, 0);
                    s = NQ(j, 0) + delh;
                    CalRho(r, s, t, rho_p);
                    s = NQ(j, 0) - delh;
                    CalRho(r, s, t, rho_m);
                    rho_s(k, all) = (rho_p - rho_m) / (2 * delh);

                    // rho_t
                    r = NQ(i, 0);
                    s = NQ(j, 0);
                    t = NQ(l, 0) + delh;
                    CalRho(r, s, t, rho_p);
                    t = NQ(l, 0) - delh;
                    CalRho(r, s, t, rho_m);
                    rho_t(k, all) = (rho_p - rho_m) / (2 * delh);

                    // PU and FE functions
                    r = NQ(i, 0);
                    s = NQ(j, 0);
                    t = NQ(l, 0);

                    CalRho(r, s, t, rho_p);
                    rho(k, all) = rho_p;

                    h(k, 0) = (1 - r) * (1 - s) * (1 - t) / 8;
                    h(k, 1) = (1 + r) * (1 - s) * (1 - t) / 8;
                    h(k, 2) = (1 + r) * (1 + s) * (1 - t) / 8;
                    h(k, 3) = (1 - r) * (1 + s) * (1 - t) / 8;
                    h(k, 4) = (1 - r) * (1 - s) * (1 + t) / 8;
                    h(k, 5) = (1 + r) * (1 - s) * (1 + t) / 8;
                    h(k, 6) = (1 + r) * (1 + s) * (1 + t) / 8;
                    h(k, 7) = (1 - r) * (1 + s) * (1 + t) / 8;

                    h_r(k, 0) = -(1 - s) * (1 - t) / 8;
                    h_r(k, 1) = (1 - s) * (1 - t) / 8;
                    h_r(k, 2) = (1 + s) * (1 - t) / 8;
                    h_r(k, 3) = -(1 + s) * (1 - t) / 8;
                    h_r(k, 4) = -(1 - s) * (1 + t) / 8;
                    h_r(k, 5) = (1 - s) * (1 + t) / 8;
                    h_r(k, 6) = (1 + s) * (1 + t) / 8;
                    h_r(k, 7) = -(1 + s) * (1 + t) / 8;

                    h_s(k, 0) = (1 - r) * (-1) * (1 - t) / 8;
                    h_s(k, 1) = (1 + r) * (-1) * (1 - t) / 8;
                    h_s(k, 2) = (1 + r) * (1 - t) / 8;
                    h_s(k, 3) = (1 - r) * (1 - t) / 8;
                    h_s(k, 4) = (1 - r) * (-1) * (1 + t) / 8;
                    h_s(k, 5) = (1 + r) * (-1) * (1 + t) / 8;
                    h_s(k, 6) = (1 + r) * (1 + t) / 8;
                    h_s(k, 7) = (1 - r) * (1 + t) / 8;

                    h_t(k, 0) = (1 - r) * (1 - s) * (-1) / 8;
                    h_t(k, 1) = (1 + r) * (1 - s) * (-1) / 8;
                    h_t(k, 2) = (1 + r) * (1 + s) * (-1) / 8;
                    h_t(k, 3) = (1 - r) * (1 + s) * (-1) / 8;
                    h_t(k, 4) = (1 - r) * (1 - s) / 8;
                    h_t(k, 5) = (1 + r) * (1 - s) / 8;
                    h_t(k, 6) = (1 + r) * (1 + s) / 8;
                    h_t(k, 7) = (1 - r) * (1 + s) / 8;
                }
            }
        }
    }

    // Upon instantiation, all basic functions are evaluated at quadrature points and saved.
    OFE8nodeLinear::OFE8nodeLinear(const std::string &type)
    {
        ////////////////////////////
        //--Function Evaluations--//
        ////////////////////////////

        if (type == "stiffness")
        {
            initialization(NQ_Stiffness);
        }
        else if (type == "mass")
        {
            initialization(NQ_Mass);
        }
        else
        {
            // Handle unexpected type
            throw std::invalid_argument("Invalid type specified");
        }
    }

    // Calculate the PU functions, Rho's
    void OFE8nodeLinear::CalRho(const double &r, const double &s, const double &t, Matrix<double, 1, 8> &rho)
    {
        rho(0) = -(0.0625) * (-1 + r) * (-1 + s) * (-1 + t) * (2 + beta * r + beta * pow(r, 2) + beta * s - 2 * beta * r * s - beta * pow(r, 2) * s + beta * pow(s, 2) - beta * r * pow(s, 2) + beta * t - 2 * beta * r * t - beta * pow(r, 2) * t - 2 * beta * s * t + 3 * beta * r * s * t + beta * pow(r, 2) * s * t - beta * pow(s, 2) * t + beta * r * pow(s, 2) * t + beta * pow(t, 2) - beta * r * pow(t, 2) - beta * s * pow(t, 2) + beta * r * s * pow(t, 2));
        rho(1) = (0.0625) * (1 + r) * (-1 + s) * (-1 + t) * (2 - beta * r + beta * pow(r, 2) + beta * s + 2 * beta * r * s - beta * pow(r, 2) * s + beta * pow(s, 2) + beta * r * pow(s, 2) + beta * t + 2 * beta * r * t - beta * pow(r, 2) * t - 2 * beta * s * t - 3 * beta * r * s * t + beta * pow(r, 2) * s * t - beta * pow(s, 2) * t - beta * r * pow(s, 2) * t + beta * pow(t, 2) + beta * r * pow(t, 2) - beta * s * pow(t, 2) - beta * r * s * pow(t, 2));
        rho(2) = (0.0625) * (1 + r) * (1 + s) * (-1 + t) * (-2 + beta * r - beta * pow(r, 2) + beta * s + 2 * beta * r * s - beta * pow(r, 2) * s - beta * pow(s, 2) - beta * r * pow(s, 2) - beta * t - 2 * beta * r * t + beta * pow(r, 2) * t - 2 * beta * s * t - 3 * beta * r * s * t + beta * pow(r, 2) * s * t + beta * pow(s, 2) * t + beta * r * pow(s, 2) * t - beta * pow(t, 2) - beta * r * pow(t, 2) - beta * s * pow(t, 2) - beta * r * s * pow(t, 2));
        rho(3) = -(0.0625) * (-1 + r) * (1 + s) * (-1 + t) * (-2 - beta * r - beta * pow(r, 2) + beta * s - 2 * beta * r * s - beta * pow(r, 2) * s - beta * pow(s, 2) + beta * r * pow(s, 2) - beta * t + 2 * beta * r * t + beta * pow(r, 2) * t - 2 * beta * s * t + 3 * beta * r * s * t + beta * pow(r, 2) * s * t + beta * pow(s, 2) * t - beta * r * pow(s, 2) * t - beta * pow(t, 2) + beta * r * pow(t, 2) - beta * s * pow(t, 2) + beta * r * s * pow(t, 2));
        rho(4) = -(0.0625) * (-1 + r) * (-1 + s) * (1 + t) * (-2 - beta * r - beta * pow(r, 2) - beta * s + 2 * beta * r * s + beta * pow(r, 2) * s - beta * pow(s, 2) + beta * r * pow(s, 2) + beta * t - 2 * beta * r * t - beta * pow(r, 2) * t - 2 * beta * s * t + 3 * beta * r * s * t + beta * pow(r, 2) * s * t - beta * pow(s, 2) * t + beta * r * pow(s, 2) * t - beta * pow(t, 2) + beta * r * pow(t, 2) + beta * s * pow(t, 2) - beta * r * s * pow(t, 2));
        rho(5) = (0.0625) * (1 + r) * (-1 + s) * (1 + t) * (-2 + beta * r - beta * pow(r, 2) - beta * s - 2 * beta * r * s + beta * pow(r, 2) * s - beta * pow(s, 2) - beta * r * pow(s, 2) + beta * t + 2 * beta * r * t - beta * pow(r, 2) * t - 2 * beta * s * t - 3 * beta * r * s * t + beta * pow(r, 2) * s * t - beta * pow(s, 2) * t - beta * r * pow(s, 2) * t - beta * pow(t, 2) - beta * r * pow(t, 2) + beta * s * pow(t, 2) + beta * r * s * pow(t, 2));
        rho(6) = (0.0625) * (1 + r) * (1 + s) * (1 + t) * (2 - beta * r + beta * pow(r, 2) - beta * s - 2 * beta * r * s + beta * pow(r, 2) * s + beta * pow(s, 2) + beta * r * pow(s, 2) - beta * t - 2 * beta * r * t + beta * pow(r, 2) * t - 2 * beta * s * t - 3 * beta * r * s * t + beta * pow(r, 2) * s * t + beta * pow(s, 2) * t + beta * r * pow(s, 2) * t + beta * pow(t, 2) + beta * r * pow(t, 2) + beta * s * pow(t, 2) + beta * r * s * pow(t, 2));
        rho(7) = -(0.0625) * (-1 + r) * (1 + s) * (1 + t) * (2 + beta * r + beta * pow(r, 2) - beta * s + 2 * beta * r * s + beta * pow(r, 2) * s + beta * pow(s, 2) - beta * r * pow(s, 2) - beta * t + 2 * beta * r * t + beta * pow(r, 2) * t - 2 * beta * s * t + 3 * beta * r * s * t + beta * pow(r, 2) * s * t + beta * pow(s, 2) * t - beta * r * pow(s, 2) * t + beta * pow(t, 2) - beta * r * pow(t, 2) + beta * s * pow(t, 2) - beta * r * s * pow(t, 2));
    }

    // OFE interpolation
    void OFE8nodeLinear::OFEintp(Matrix<double, 1, 32> &H, const double &x, const double &y, const double &z, const int &k, const MatrixXd &nodes, const Matrix<int, 1, 8> &element, const std::vector<double> &radius)
    {
        for (int j = 0; j < 8; j++)
        {
            xK = (x - nodes(element[j], 0)) / radius[element[j]];
            yK = (y - nodes(element[j], 1)) / radius[element[j]];
            zK = (z - nodes(element[j], 2)) / radius[element[j]];
            H(4 * j) = rho(k, j);
            H(4 * j + 1) = rho(k, j) * xK;
            H(4 * j + 2) = rho(k, j) * yK;
            H(4 * j + 3) = rho(k, j) * zK;
        }
    }
    void OFE8nodeLinear::ElementMass(const double &density, const MatrixXd &nodes, const Matrix<int, 1, 8> &element, std::vector<Triplet<double>> &triplets, const std::vector<double> &radius)
    {
        M.setZero();
        int nodeNum = nodes.rows();

        for (int i = 0, k = 0; i < NQ_Mass.rows(); i++)
        {
            for (int j = 0; j < NQ_Mass.rows(); j++)
            {
                for (int l = 0; l < NQ_Mass.rows(); l++, k++)
                {
                    Jacobian(nodes, element, k);
                    detJ = J.determinant();
                    x = h(k, all).dot(nodes(element, 0));
                    y = h(k, all).dot(nodes(element, 1));
                    z = h(k, all).dot(nodes(element, 2));
                    OFEintp(H, x, y, z, k, nodes, element, radius);
                    M.noalias() += H.transpose() * H * detJ * weight(k) * density;
                }
            }
        }

        std::vector<int> idx;
        for (int i = 0; i < 8; i++)
        {
            idx.push_back(4 * (element[i]));
            idx.push_back(4 * (element[i]) + 1);
            idx.push_back(4 * (element[i]) + 2);
            idx.push_back(4 * (element[i]) + 3);
        }

        for (int i = 0; i < M.rows(); i++)
        {
            for (int j = 0; j < M.cols(); j++)
            {
                Triplet<double> tr_x(idx[i] * 3, idx[j] * 3, M(i, j));
                triplets.push_back(tr_x);
                Triplet<double> tr_y(idx[i] * 3 + 1, idx[j] * 3 + 1, M(i, j));
                triplets.push_back(tr_y);
                Triplet<double> tr_z(idx[i] * 3 + 2, idx[j] * 3 + 2, M(i, j));
                triplets.push_back(tr_z);
            }
        }
    }

    void OFE8nodeLinear::ElementStiffness(const MatrixXd &C, const MatrixXd &nodes, const Matrix<int, 1, 8> &element, std::vector<Triplet<double>> &triplets, const std::vector<double> &radius) //, const std::vector<double>& radius
    {

        // iteration over quadrature points
        K.setZero();
        int nodeNum = nodes.rows();
        for (int i = 0, k = 0; i < NQ_Stiffness.rows(); i++)
        {
            for (int s = 0; s < NQ_Stiffness.rows(); s++)
            {
                for (int l = 0; l < NQ_Stiffness.rows(); l++, k++)
                {
                    Jacobian(nodes, element, k);
                    invJ = J.inverse();
                    detJ = J.determinant();

                    x = h(k, all).dot(nodes(element, 0));
                    y = h(k, all).dot(nodes(element, 1));
                    z = h(k, all).dot(nodes(element, 2));
                    // iteration over nodes to obtain dHdx, dHdy, and dHdz
                    for (int j = 0; j < 8; j++)
                    {
                        xK = (x - nodes(element[j], 0)) / radius[element[j]];
                        yK = (y - nodes(element[j], 1)) / radius[element[j]];
                        zK = (z - nodes(element[j], 2)) / radius[element[j]];
                        rhoK_x = rho_r(k, j) * invJ(0, 0) + rho_s(k, j) * invJ(0, 1) + rho_t(k, j) * invJ(0, 2);
                        rhoK_y = rho_r(k, j) * invJ(1, 0) + rho_s(k, j) * invJ(1, 1) + rho_t(k, j) * invJ(1, 2);
                        rhoK_z = rho_r(k, j) * invJ(2, 0) + rho_s(k, j) * invJ(2, 1) + rho_t(k, j) * invJ(2, 2);

                        dHdx(4 * j) = rhoK_x;
                        dHdx(4 * j + 1) = rhoK_x * xK + rho(k, j) / radius[element[j]];
                        dHdx(4 * j + 2) = rhoK_x * yK;
                        dHdx(4 * j + 3) = rhoK_x * zK;

                        dHdy(4 * j) = rhoK_y;
                        dHdy(4 * j + 1) = rhoK_y * xK;
                        dHdy(4 * j + 2) = rhoK_y * yK + rho(k, j) / radius[element[j]];
                        dHdy(4 * j + 3) = rhoK_y * zK;

                        dHdz(4 * j) = rhoK_z;
                        dHdz(4 * j + 1) = rhoK_z * xK;
                        dHdz(4 * j + 2) = rhoK_z * yK;
                        dHdz(4 * j + 3) = rhoK_z * zK + rho(k, j) / radius[element[j]];
                    }
                    B(0, seq(0, 31)) = dHdx;
                    B(1, seq(32, 63)) = dHdy;
                    B(2, seq(64, 95)) = dHdz;

                    B(3, seq(0, 31)) = dHdy;
                    B(3, seq(32, 63)) = dHdx;

                    B(4, seq(32, 63)) = dHdz;
                    B(4, seq(64, 95)) = dHdy;

                    B(5, seq(0, 31)) = dHdz;
                    B(5, seq(64, 95)) = dHdx;

                    K.noalias() += B.transpose() * C * B * detJ * weight(k);
                }
            }
        }
        // std::cout << K << std::endl;

        std::vector<int> idx;
        for (int i = 0; i < 8; i++)
        {
            idx.push_back(4 * (element[i]));
            idx.push_back(4 * (element[i]) + 1);
            idx.push_back(4 * (element[i]) + 2);
            idx.push_back(4 * (element[i]) + 3);
        }

        // int idx[9];
        // for (int i = 0; i < 3; i++)
        // {
        //     idx[3 * i] = 3 * (element[i]);
        //     idx[3 * i + 1] = 3 * (element[i]) + 1;
        //     idx[3 * i + 2] = 3 * (element[i]) + 2;
        // }

        for (int i = 0; i < K.rows(); i++)
        {
            for (int j = 0; j < K.cols(); j++)
            {
                if (i < 32)
                    row = 3 * idx[i]; // idx[i]; // element[i];
                else if (i > 31 & i < 64)
                    row = 3 * idx[i - 32] + 1; // idx[i - 32] + 4 * nodeNum;
                else if (i > 63 & i < 96)
                    row = 3 * idx[i - 64] + 2; // idx[i - 64] + 8 * nodeNum;
                else
                    std::cout
                        << "Error occurs with index assignment (row)" << std::endl;

                if (j < 32)
                    col = 3 * idx[j]; // idx[j];
                else if (j > 31 & j < 64)
                    col = 3 * idx[j - 32] + 1; // idx[j - 32] + 4 * nodeNum;
                else if (j > 63 & j < 96)
                    col = 3 * idx[j - 64] + 2; // idx[j - 64] + 8 * nodeNum;
                else
                    std::cout << "Error occurs with index assignment (col)" << std::endl;

                Triplet<double> tr(row, col, K(i, j));

                triplets.push_back(tr);
            }
        }
    }

    void OFE8nodeLinear::Jacobian(const MatrixXd &nodes, const Matrix<int, 1, 8> &element, const int &i)
    {
        J << h_r(i, all) * nodes(element, 0), h_r(i, all) * nodes(element, 1), h_r(i, all) * nodes(element, 2),
            h_s(i, all) * nodes(element, 0), h_s(i, all) * nodes(element, 1), h_s(i, all) * nodes(element, 2),
            h_t(i, all) * nodes(element, 0), h_t(i, all) * nodes(element, 1), h_t(i, all) * nodes(element, 2);
    }

    std::vector<double> OFE8nodeLinear::radius(const MatrixXd &nodes, const Eigen::MatrixXi &elements, const size_t &ne, const size_t &nn)
    {

        // find the combinations
        int arr[] = {0, 1, 2, 3, 4, 5, 6, 7};
        std::vector<std::array<int, 2>> combi;
        Combination CalCombi;
        CalCombi.getCombination(arr, 8, 2, combi);

        std::vector<std::vector<double>> radii(nn);
        double X0, X1, Y0, Y1, Z0, Z1, d;
        for (int i = 0; i < ne; i++)
        {

            for (auto it = begin(combi); it != end(combi); ++it)
            {
                X0 = nodes(elements(i, it->at(0)), 0);
                Y0 = nodes(elements(i, it->at(0)), 1);
                Z0 = nodes(elements(i, it->at(0)), 2);
                X1 = nodes(elements(i, it->at(1)), 0);
                Y1 = nodes(elements(i, it->at(1)), 1);
                Z1 = nodes(elements(i, it->at(1)), 2);
                d = sqrt(pow(X0 - X1, 2.0) + pow(Y0 - Y1, 2.0) + pow(Z0 - Z1, 2.0));
                radii[elements(i, it->at(0))].push_back(d);
                radii[elements(i, it->at(1))].push_back(d);
            }

            // double X0 = nodes(elements[i].nodesIds[0], 0);
            // double Y0 = nodes(elements[i].nodesIds[0], 1);
            // double X1 = nodes(elements[i].nodesIds[1], 0);
            // double Y1 = nodes(elements[i].nodesIds[1], 1);
            // double X2 = nodes(elements[i].nodesIds[2], 0);
            // double Y2 = nodes(elements[i].nodesIds[2], 1);

            // double d01 = sqrt(pow(X0 - X1, 2.0) + pow(Y0 - Y1, 2.0));
            // double d12 = sqrt(pow(X2 - X1, 2.0) + pow(Y2 - Y1, 2.0));
            // double d20 = sqrt(pow(X2 - X0, 2.0) + pow(Y2 - Y0, 2.0));
            // radii[elements[i].nodesIds[0]].push_back(d01);
            // radii[elements[i].nodesIds[0]].push_back(d20);
            // radii[elements[i].nodesIds[1]].push_back(d01);
            // radii[elements[i].nodesIds[1]].push_back(d12);
            // radii[elements[i].nodesIds[2]].push_back(d20);
            // radii[elements[i].nodesIds[2]].push_back(d12);
        }
        std::vector<double> radius;
        for (auto iter = radii.begin(); iter != radii.end(); iter++)
        {
            radius.push_back(*max_element(iter->begin(), iter->end()) / 2);
        }
        return radius;
    }
}