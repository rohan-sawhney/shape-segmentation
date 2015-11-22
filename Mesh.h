#ifndef MESH_H
#define MESH_H

#include "Types.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"

class Mesh {
public:
    // default constructor
    Mesh();
        
    // read mesh from file
    bool read(const std::string& fileName);
    
    // write mesh to file
    bool write(const std::string& fileName) const;
    
    // segments shape using the shape diameter function
    void segment(const int clusters);
            
    // member variables
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> vertices;
    std::vector<Eigen::Vector3d> uvs;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<HalfEdgeIter> boundaries;

private:
    // center mesh about origin and rescale to unit radius
    void normalize();
    
    // computes normalized sdf for each face
    void computeNormalizedSDFValues();
};

#endif