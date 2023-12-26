#include <stdexcept>

#include "Nodes.h"
#include "PreProcessing.h"
namespace Dynamis::core
{
    Nodes::Nodes(std::unique_ptr<std::ifstream> &file)
    {
        if (!Dynamis::PreProcessing::readValue(file, "NumNodes", numNodes))
        {
            throw std::runtime_error("Error in reading NumNodes");
        }
        else
        {
            std::cout << "Number of nodes: " << numNodes << std::endl;
        }

        if (!Dynamis::PreProcessing::readMatrix(file, "Nodes", nodes, numNodes, 3))
        {
            throw std::runtime_error("Error in reading nodes");
        }
    }
    Eigen::MatrixXd &Nodes::getNodes()
    {
        return nodes;
    }
    size_t &Nodes::getNumNodes()
    {
        return numNodes;
    }
}