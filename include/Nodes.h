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
        size_t &getNumNodes();

    private:
        // Member variables
        Eigen::MatrixXd nodes;
        size_t numNodes;
    };

} // namespace Dynamis::core

#endif // NODES_H
