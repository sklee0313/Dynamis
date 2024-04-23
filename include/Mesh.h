#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <fstream>
#include <memory>
#include "Eigen/Dense"

namespace Dynamis::core
{

    class Mesh
    {
    public:
        // Constructor
        Mesh(std::unique_ptr<std::ifstream> &file);

        // Member functions
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &getNodes();
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &getElements();
        size_t &getNumNodes();
        size_t &getNumElements();
        size_t &getNumVertices();

    private:
        // Member variables
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes;
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> elements;
        size_t numNodes;
        size_t numVertices;
        size_t numElements;
    };

} // namespace Dynamis::core

#endif // MESH_H
