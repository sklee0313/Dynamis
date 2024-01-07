#include <stdexcept>

#include "Mesh.h"
#include "PreProcessing.h"
namespace Dynamis::core
{
    Mesh::Mesh(std::unique_ptr<std::ifstream> &file)
    {
        ///////////////////////
        //--Read nodes data--//
        ///////////////////////

        if (!Dynamis::PreProcessing::readValue(file, "NumNodes", numNodes))
        {
            throw std::runtime_error("Error in reading NumNodes");
        }
        else
        {
            std::cout << "Number of nodes: " << numNodes << std::endl;
        }

        if (!Dynamis::PreProcessing::readMatrix<double>(file, "Nodes", nodes, numNodes, 3))
        {
            throw std::runtime_error("Error in reading nodes");
        }

        //////////////////////////
        //--Read elements data--//
        //////////////////////////

        if (!Dynamis::PreProcessing::readValue(file, "NumElements", numElements))
        {
            throw std::runtime_error("Error in reading NumElements");
        }
        else
        {
            std::cout << "Number of elements: " << numElements << std::endl;
        }

        if (!Dynamis::PreProcessing::readValue(file, "NumVertices", numVertices))
        {
            throw std::runtime_error("Error in reading NumEVertices");
        }
        else
        {
            std::cout << "Number of vertices: " << numVertices << std::endl;
        }

        if (!Dynamis::PreProcessing::readMatrix<int>(file, "Elements", elements, numElements, numVertices))
        {
            throw std::runtime_error("Error in reading elements");
        }
    }
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &Mesh::getNodes()
    {
        return nodes;
    }
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &Mesh::getElements()
    {
        return elements;
    }
    size_t &Mesh::getNumNodes()
    {
        return numNodes;
    }
    size_t &Mesh::getNumElements()
    {
        return numElements;
    }
    size_t &Mesh::getNumVertices()
    {
        return numVertices;
    }
}