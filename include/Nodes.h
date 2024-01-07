#ifndef NODES_H
#define NODES_H

#include <iostream>
#include <fstream>
#include <memory>
#include "Eigen/Dense"

namespace Dynamis::core
{

    class Nodes
    {
    public:
        // Constructor
        Nodes(std::unique_ptr<std::ifstream> &file);

        // Member functions
        Eigen::MatrixXd &getNodes();
        Eigen::MatrixXd &getElements();
        size_t &getNumNodes();
        size_t &getNumElements();
        size_t &getNumVertices();

    private:
        // Member variables
        Eigen::MatrixXd nodes;
        Eigen::MatrixXd elements;
        size_t numNodes;
        size_t numVertices;
        size_t numElements;
    };

} // namespace Dynamis::core

#endif // NODES_H
