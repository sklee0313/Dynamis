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
        std::vector<Node> nodes; // set of nodes
        size_t nn;               // number of nodes

    public:
        Nodes(std::ifstream &file);
    };
}

#endif // NODES_H
