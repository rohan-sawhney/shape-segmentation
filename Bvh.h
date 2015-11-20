#ifndef BVH_H
#define BVH_H

#include "Types.h"
#include "BoundingBox.h"

class Node {
public:
    // member variables
    BoundingBox boundingBox;
    int startId, range, rightOffset;
};

class Bvh {
public:
    Bvh(Mesh *meshPtr0, const int leafSize0 = 1);
    
    // returns distance to closest intersection point
    double distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& normal) const;
    
private:    
    
    // builds the bvh
    void build();
    
    int nodeCount, leafCount, leafSize;
    std::vector<Node> flatTree;
    Mesh *meshPtr;
};

#endif
