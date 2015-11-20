#include "Bvh.h"
#include "Mesh.h"
#include <stack>

Bvh::Bvh(Mesh *meshPtr0, const int leafSize0):
meshPtr(meshPtr0),
nodeCount(0),
leafCount(0),
leafSize(leafSize0)
{
    build();
}

struct NodeEntry {
    // parent id
    int parentId;
    
    // range of objects covered by the node
    int startId, endId;
};

class TraversalEntry {
public:
    // constructor
    TraversalEntry(const int id0, const double d0): id(id0), d(d0) {}
    
    // id
    int id;
    
    // distance
    double d;
};

void Bvh::build()
{
    std::stack<NodeEntry> stack;
    const int faceCount = (int)meshPtr->faces.size();
    
    NodeEntry nodeEntry;
    nodeEntry.parentId = -1;
    nodeEntry.startId = 0;
    nodeEntry.endId = faceCount;
    stack.push(nodeEntry);
    
    Node node;
    std::vector<Node> nodes;
    nodes.reserve(faceCount * 2);
    
    while (!stack.empty()) {
        // pop item off the stack and create a node
        nodeEntry = stack.top();
        stack.pop();
        int startId = nodeEntry.startId;
        int endId = nodeEntry.endId;
        
        nodeCount ++;
        node.startId = startId;
        node.range = endId - startId;
        node.rightOffset = 2;
        
        // calculate bounding box
        BoundingBox boundingBox(meshPtr->faces[startId].boundingBox());
        BoundingBox boundingCentroid(meshPtr->faces[startId].centroid());
        for (int i = startId+1; i < endId; i++) {
            boundingBox.expandToInclude(meshPtr->faces[i].boundingBox());
            boundingCentroid.expandToInclude(meshPtr->faces[i].centroid());
        }
        node.boundingBox = boundingBox;
        
        // if node is a leaf
        if (node.range <= leafSize) {
            node.rightOffset = 0;
            leafCount ++;
        }
        
        nodes.push_back(node);
        
        // compute parent's rightOffset
        if (nodeEntry.parentId != -1) {
            nodes[nodeEntry.parentId].rightOffset --;
            
            if (nodes[nodeEntry.parentId].rightOffset == 0) {
                nodes[nodeEntry.parentId].rightOffset = nodeCount - 1 - nodeEntry.parentId;
            }
        }
        
        // if a leaf, no need to subdivide
        if (node.rightOffset == 0) {
            continue;
        }
        
        // find the center of the longest dimension
        int maxDimension = boundingCentroid.maxDimension();
        double splitCoord = 0.5 * (boundingCentroid.min[maxDimension] +
                                   boundingCentroid.max[maxDimension]);
        
        // partition faces
        int mid = startId;
        for (int i = startId; i < endId; i++) {
            if (meshPtr->faces[i].centroid()[maxDimension] < splitCoord) {
                std::swap(meshPtr->faces[i], meshPtr->faces[mid]);
                mid ++;
            }
        }
        
        // in case of a bad split
        if (mid == startId || mid == endId) {
            mid = startId + (endId-startId) / 2;
        }
        
        // push right child
        nodeEntry.startId = mid;
        nodeEntry.endId = endId;
        nodeEntry.parentId = nodeCount - 1;
        stack.push(nodeEntry);
        
        // push left child
        nodeEntry.startId = startId;
        nodeEntry.endId = mid;
        nodeEntry.parentId = nodeCount - 1;
        stack.push(nodeEntry);
    }
    
    // copy node data into temp array
    for (int i = 0; i < nodeCount; i ++) {
        flatTree.push_back(nodes[i]);
    }
}

double Bvh::distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& normal) const
{
    double minD = INFINITY;
    
    int id = 0;
    int closer, further;
    double dist1 = 0.0;
    double dist2 = 0.0;
    
    TraversalEntry t(id, -INFINITY);
    std::stack<TraversalEntry> stack;
    stack.push(t);
    
    while (!stack.empty()) {
        TraversalEntry t = stack.top();
        id = t.id;
        stack.pop();
        
        if (minD < t.d) continue;
        
        const Node &node(flatTree[id]);
        // node is a leaf
        if (node.rightOffset == 0) {
            for (int i = 0; i < node.range; i++) {
                // check for overlap
                double d = meshPtr->faces[node.startId+i].distance(origin, normal);
                if (d < minD) {
                    minD = d;
                }
            }
            
        } else { // not a leaf
            bool hit0 = flatTree[id+1].boundingBox.intersect(origin, normal, dist1);
            bool hit1 = flatTree[id+node.rightOffset].boundingBox.intersect(origin, normal, dist2);
            
            // hit both bounding boxes
            if (hit0 && hit1) {
                closer = id+1;
                further = id+node.rightOffset;
                
                if (dist2 < dist1) {
                    std::swap(dist1, dist2);
                    std::swap(closer, further);
                }
                
                // push farther node first
                stack.push(TraversalEntry(further, dist2));
                stack.push(TraversalEntry(closer, dist1));
                
            } else if (hit0) {
                stack.push(TraversalEntry(id+1, dist1));
                
            } else if (hit1) {
                stack.push(TraversalEntry(id+node.rightOffset, dist2));
            }
        }
    }

    return minD;
}
