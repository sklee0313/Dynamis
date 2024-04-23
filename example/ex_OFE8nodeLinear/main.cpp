#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <stack>
#include <ctime>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Spectra/MatOp/SparseGenMatProd.h"
#include "Spectra/SymGEigsShiftSolver.h"

#include "igl/slice.h"
#include "igl/setdiff.h"
#include "OFE8nodeLinear.h"
#include "Element.h"
#include "RCMordering.h"
#include "TicToc.h"
#include "PreProcessing.h"
#include "LinearElasticity.h"
#include "Mesh.h"

using namespace Eigen;

int main(int argc, char *argv[])
{
    ///////////////////////////////
    //--Check if input is valid--//
    ///////////////////////////////
    std::cout << "hello " << std::endl;

    auto infile = Dynamis::PreProcessing::InputChecker(argc, argv); // fstream returned

    ////////////////////////////
    //--3D Linear Elasticity--//
    ////////////////////////////

    Dynamis::ConstitutiveLaw::LinearElasticity elasticity(infile); // read material data and instantiate the law

    const auto C = elasticity.getC();
    const auto density = elasticity.getDensity();

    Dynamis::core::Mesh mesh(infile);
    const Eigen::MatrixXd &nodes = mesh.getNodes();
    const size_t &nn = mesh.getNumNodes();
    const size_t &ne = mesh.getNumElements();
    const size_t &npe = mesh.getNumVertices();
    const MatrixXi &elements = mesh.getElements();

    //////////////////////
    //--Discretization--//
    //////////////////////

    Dynamis::OFE8nodeLinear ElementMethod_K = Dynamis::OFE8nodeLinear("stiffness");
    Dynamis::OFE8nodeLinear ElementMethod_M = Dynamis::OFE8nodeLinear("mass");

    std::vector<double> radius = ElementMethod_K.radius(nodes, elements, ne, nn);
    // Construction of Stiffness Matrix using Direct Stiffness Method
    tic();
    std::vector<Triplet<double>> triplets;
    for (int i = 0; i < ne; i++)
    {
        ElementMethod_K.ElementStiffness(C, nodes, elements(i, all), triplets, radius); //,NQ
    }

    SparseMatrix<double> globalK(3 * 4 * nn, 3 * 4 * nn);
    globalK.setFromTriplets(triplets.begin(), triplets.end());
    toc();
    std::cout << "for stiffness matrix construction" << std::endl;

    // Construction of Mass Matrix
    tic();
    triplets.clear();
    for (int i = 0; i < ne; i++)
    {
        ElementMethod_M.ElementMass(density, nodes, elements(i, all), triplets, radius); //,NQ
    }

    SparseMatrix<double> globalM(3 * 4 * nn, 3 * 4 * nn);
    globalM.setFromTriplets(triplets.begin(), triplets.end());
    toc();
    std::cout << "for mass matrix construction" << std::endl;

    ///////////////////////////////
    // apply boundary conditions //
    ///////////////////////////////

    std::vector<int> II;
    for (int i = 0; i < nn; i++)
    {
        if ((abs(nodes(i, 1)) < 1e-5 & nodes(i, 0) < 0)) // if the nodal y-position is -1
        {
            II.push_back(ElementMethod_M.nodaldofs * i); // we kill 1, x, z
            II.push_back(ElementMethod_M.nodaldofs * i + 1);
            II.push_back(ElementMethod_M.nodaldofs * i + 2);

            II.push_back(ElementMethod_M.nodaldofs * i + 3);
            II.push_back(ElementMethod_M.nodaldofs * i + 4);
            II.push_back(ElementMethod_M.nodaldofs * i + 5);

            II.push_back(ElementMethod_M.nodaldofs * i + 9);
            II.push_back(ElementMethod_M.nodaldofs * i + 10);
            II.push_back(ElementMethod_M.nodaldofs * i + 11);
        }
    }
    VectorXi FixedDofsIndices = VectorXi::Map(II.data(), static_cast<int>(II.size()));
    VectorXi AllIndices = VectorXi::LinSpaced(ElementMethod_M.nodaldofs * nn, 0, ElementMethod_M.nodaldofs * nn);
    VectorXi FreeDofsIndices, IA;
    igl::setdiff(AllIndices, FixedDofsIndices, FreeDofsIndices, IA);

    SparseMatrix<double> globalKK(FreeDofsIndices.size(), FreeDofsIndices.size());
    SparseMatrix<double> globalMM(FreeDofsIndices.size(), FreeDofsIndices.size());
    igl::slice(globalK, FreeDofsIndices, FreeDofsIndices, globalKK);
    igl::slice(globalM, FreeDofsIndices, FreeDofsIndices, globalMM);
    std::cout << "dofs before imposing boundary condtions = " << globalK.rows() << std::endl;
    std::cout << "dofs after imposing boundary condtions = " << globalKK.rows() << std::endl;

    ReorderingSSM m(globalKK);
    VectorXi r = m.ReverseCuthillMckee();
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
    perm.indices() = r;
    globalKK = perm.transpose() * globalKK * perm;
    globalMM = perm.transpose() * globalMM * perm;

    ////////////////////////////////////
    // generalized eigenvalue problem //
    ////////////////////////////////////

    // The generalized eigenvalue problem is solving using Spectra

    // Construct matrix operation objects using the wrapper classes
    using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = Spectra::SparseSymMatProd<double>;
    OpType op(globalKK, globalMM);
    BOpType Bop(globalMM);

    // seeking 10 eigenvalues
    Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert>
        geigs(op, Bop, 10, 20, 0.0);

    // Initialize and compute
    tic();
    geigs.init();
    int nconv = geigs.compute(Spectra::SortRule::LargestMagn);
    std::cout << "Solution of the generalized eigenvalue problem takes ";
    toc();

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;
    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Number of converged generalized eigenvalues: " << nconv << std::endl;
    std::cout << "Natural frequencies found (Hz) :\n"
              << evalues.reverse().cwiseSqrt() / (2 * 3.141592) << std::endl;

    return 0; // Return 0 to indicate successful execution.
}
