#ifndef NODE_H
#define NODE_H

#include <Eigen/Dense>

class Node
{
private:
    Eigen::Vector3d position; // Position vector (x, y, z)

public:
    // Constructor to initialize the position
    Node(double x, double y, double z) : position(Eigen::Vector3d(x, y, z)) {}

    // Accessor methods to get the position data
    double getX() const { return position.x(); }
    double getY() const { return position.y(); }
    double getZ() const { return position.z(); }

    // Additional methods as needed
};

#endif // NODE_H