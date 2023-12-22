#include <stdexcept>

#include "Nodes.h"
#include "PreProcessing.h"
namespace Dynamis::core
{
    Nodes::Nodes(std::ifstream &file)
    {
        bool numNodesRead = false, nodesRead = false;

        std::string line;
        file.seekg(0, std::ios::beg); // Reset to the beginning of the file
        while (std::getline(file, line) && !numNodesRead)
        {
            if (!numNodesRead && Dynamis::PreProcessing::tryReadValue(file, line, "NumNodes", numNodes))
            {
                numNodesRead = true;
                std::cout << "Number of nodes: " << numNodes << std::endl;
                break;
            }
        }
        if (!numNodesRead)
        {
            throw std::runtime_error("Error in reading number of nodes");
        }
        file.seekg(0, std::ios::beg); // Reset to the beginning of the file
        while (std::getline(file, line) && !nodesRead)
        {
            if (!nodesRead && Dynamis::PreProcessing::tryReadMatrix(file, line, "Nodes", nodes, numNodes, 3))
            {
                numNodesRead = true;
                std::cout << "Nodes creation done" << std::endl;
                break;
            }
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