#ifndef HALF_EDGE_H
#define HALF_EDGE_H

#include "Types.h"

class HalfEdge {
public:
    // next halfedge around the current face
    HalfEdgeIter next;
    
    // other halfedge associated with this edge
    HalfEdgeIter flip;
    
    // vertex at the tail of the halfedge
    VertexIter vertex;
    
    // edge associated with this halfedge
    EdgeIter edge;
    
    // face associated with this halfedge
    FaceIter face;
    
    // id between 0 and |HE|-1
    int index;
    
    // uv associated with vertex at tail of halfedge
    Eigen::Vector3d uv;
    
    // normal associated with vertex at tail of halfedge
    Eigen::Vector3d normal;
    
    // checks if this halfedge is contained in boundary loop
    bool onBoundary;
};

#endif
