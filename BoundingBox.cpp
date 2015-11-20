#include "BoundingBox.h"
#define EPSILON 1e-6

BoundingBox::BoundingBox():
min(Eigen::Vector3d::Zero()),
max(Eigen::Vector3d::Zero()),
extent(Eigen::Vector3d::Zero())
{
    
}

BoundingBox::BoundingBox(const Eigen::Vector3d& min0, const Eigen::Vector3d& max0):
min(min0),
max(max0)
{
    extent = max - min;
}

BoundingBox::BoundingBox(const Eigen::Vector3d& p):
min(p),
max(p)
{
    extent = max - min;
}

void BoundingBox::expandToInclude(const Eigen::Vector3d& p)
{
    if (min.x() > p.x()) min.x() = p.x();
    if (min.y() > p.y()) min.y() = p.y();
    if (min.z() > p.z()) min.z() = p.z();
    
    if (max.x() < p.x()) max.x() = p.x();
    if (max.y() < p.y()) max.y() = p.y();
    if (max.z() < p.z()) max.z() = p.z();
    
    extent = max - min;
}

void BoundingBox::expandToInclude(const BoundingBox& b)
{
    if (min.x() > b.min.x()) min.x() = b.min.x();
    if (min.y() > b.min.y()) min.y() = b.min.y();
    if (min.z() > b.min.z()) min.z() = b.min.z();
    
    if (max.x() < b.max.x()) max.x() = b.max.x();
    if (max.y() < b.max.y()) max.y() = b.max.y();
    if (max.z() < b.max.z()) max.z() = b.max.z();
    
    extent = max - min;
}

int BoundingBox::maxDimension() const
{
    int result = 0;
    if (extent.y() > extent.x()) result = 1;
    if (extent.z() > extent.y() && extent.z() > extent.x()) result = 2;
    
    return result;
}

bool BoundingBox::intersect(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction,
                            double& dist) const
{
    double ox = origin.x();
    double dx = direction.x();
    double tMin, tMax;
    if (dx >= 0) {
        tMin = (min.x() - ox) / dx;
        tMax = (max.x() - ox) / dx;
        
    } else {
        tMin = (max.x() - ox) / dx;
        tMax = (min.x() - ox) / dx;
    }
    
    double oy = origin.y();
    double dy = direction.y();
    double tyMin, tyMax;
    if (dy >= 0) {
        tyMin = (min.y() - oy) / dy;
        tyMax = (max.y() - oy) / dy;
        
    } else {
        tyMin = (max.y() - oy) / dy;
        tyMax = (min.y() - oy) / dy;
    }
    
    if (tMin > tyMax || tyMin > tMax) {
        return false;
    }
    
    if (tyMin > tMin) tMin = tyMin;
    if (tyMax < tMax) tMax = tyMax;
    
    double oz = origin.z();
    double dz = direction.z();
    double tzMin, tzMax;
    if (dz >= 0) {
        tzMin = (min.z() - oz) / dz;
        tzMax = (max.z() - oz) / dz;
        
    } else {
        tzMin = (max.z() - oz) / dz;
        tzMax = (min.z() - oz) / dz;
    }
    
    if (tMin > tzMax || tzMin > tMax) {
        return false;
    }
    
    if (tzMin > tMin) tMin = tzMin;
    if (tzMax < tMax) tMax = tzMax;
    
    dist = tMin;
    return true;
}
