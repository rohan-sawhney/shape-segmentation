#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"
#include "BoundingBox.h"
#define EPSILON 1e-6

bool Face::isBoundary() const
{
    return he->onBoundary;
}

double Face::area() const
{
    if (isBoundary()) {
        return 0;
    }
    
    return 0.5 * normal().norm();
}

Eigen::Vector3d Face::normal() const
{
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    Eigen::Vector3d v1 = b - a;
    Eigen::Vector3d v2 = c - a;
    
    return v1.cross(v2);
}

BoundingBox Face::boundingBox() const
{
    if (isBoundary()) {
        return BoundingBox(he->vertex->position, he->next->vertex->position);
    }
    
    const Eigen::Vector3d& p1(he->vertex->position);
    const Eigen::Vector3d& p2(he->next->vertex->position);
    const Eigen::Vector3d& p3(he->next->next->vertex->position);
    
    Eigen::Vector3d min = p1;
    Eigen::Vector3d max = p1;
    
    if (p2.x() < min.x()) min.x() = p2.x();
    if (p3.x() < min.x()) min.x() = p3.x();
    
    if (p2.y() < min.y()) min.y() = p2.y();
    if (p3.y() < min.y()) min.y() = p3.y();
    
    if (p2.z() < min.z()) min.z() = p2.z();
    if (p3.z() < min.z()) min.z() = p3.z();
    
    if (p2.x() > max.x()) max.x() = p2.x();
    if (p3.x() > max.x()) max.x() = p3.x();
    
    if (p2.y() > max.y()) max.y() = p2.y();
    if (p3.y() > max.y()) max.y() = p3.y();
    
    if (p2.z() > max.z()) max.z() = p2.z();
    if (p3.z() > max.z()) max.z() = p3.z();
    
    return BoundingBox(min, max);
}

Eigen::Vector3d Face::centroid() const
{
    Eigen::Vector3d centroid;
    
    if (isBoundary()) {
        centroid = (he->vertex->position +
                    he->next->vertex->position) / 2.0;
        
    } else {
        centroid = (he->vertex->position +
                    he->next->vertex->position +
                    he->next->next->vertex->position) / 3.0;
    }
    
    return centroid;
}

bool sameDirection(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
{
    double c = v1.dot(v2);
    if (c < -1.0) c = -1.0;
    else if (c >  1.0) c = 1.0;
    
    double angle = acos(c) * 180.0 / M_PI;
    if (angle <= 90.0) {
        return true;
    }
    
    return false;
}

double Face::distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction) const
{
    // check for false intersection with the outside of the mesh
    if (!sameDirection(direction, normal().normalized())) {
        return INFINITY;
    }

    // Möller–Trumbore intersection algorithm
    const Eigen::Vector3d& p1(he->vertex->position);
    const Eigen::Vector3d& p2(he->next->vertex->position);
    const Eigen::Vector3d& p3(he->next->next->vertex->position);
    
    Eigen::Vector3d e1 = p2 - p1;
    Eigen::Vector3d e2 = p3 - p1;
    Eigen::Vector3d n = direction.cross(e2);
    
    double det = e1.dot(n);
    // ray does not lie in the plane
    if (det > -EPSILON && det < EPSILON) {
        return INFINITY;
    }
    
    double invDet = 1.0 / det;
    Eigen::Vector3d t = origin - p1;
    double u = t.dot(n) * invDet;
    
    // ray lies outside triangle
    if (u < 0.0 || u > 1.0) {
        return INFINITY;
    }
    
    Eigen::Vector3d q = t.cross(e1);
    double v = direction.dot(q) * invDet;
    // ray lies outside the triangle
    if (v < 0.0 || v + u > 1.0) {
        return INFINITY;
    }
    
    double s = e2.dot(q) * invDet;
    // intersection
    if (s > EPSILON) {
        return s;
    }
    
    // no hit
    return INFINITY;
}
