#ifndef NODES_H
#define NODES_H

#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

#include "Node.h"

namespace Dynamis::core
{

    class Nodes
    {
    private:
        Eigen::MatrixXd nodes; // set of nodes
        size_t numNodes;       // number of nodes

    public:
        Nodes(std::ifstream &file);
        Eigen::MatrixXd &getNodes();
        size_t &getNumNodes();
    };
}

#endif // NODES_H
