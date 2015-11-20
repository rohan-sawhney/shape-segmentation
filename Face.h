#ifndef FACE_H
#define FACE_H

#include "Types.h"
class BoundingBox;

class Face {
public:
    // one of the halfedges associated with this face
    HalfEdgeIter he;
    
    // id between 0 and |F|-1
    int index;
    
    // sdf va
    double sdf;
    
    // cluster
    int cluster;
        
    // checks if this face lies on boundary
    bool isBoundary() const;
    
    // returns face area
    double area() const;
    
    // computes the bounding box of the face
    BoundingBox boundingBox() const;
    
    // computes the centroid of the face
    Eigen::Vector3d centroid() const;
    
    // returns distance from origin to face
    double distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction) const;
    
    // returns normal to face
    Eigen::Vector3d normal() const;
};

#endif
