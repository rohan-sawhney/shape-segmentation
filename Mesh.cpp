#include "Mesh.h"
#include "MeshIO.h"
#include "Bvh.h"
#include <algorithm>

Mesh::Mesh()
{
    
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        normalize();
        computeNormalizedSDFValues();
    }
    
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    return false;
}

double randomDirectionInCone(Eigen::Vector3d& newDirection, const Eigen::Vector3d& direction)
{
    static const double maxAngle = 60.0 * M_PI / 180.0;
    double r = (double)rand() / RAND_MAX;
    double angle = -maxAngle + r * 2 * maxAngle;
    double l = tan(angle);
    
    Eigen::Vector3d randVec = Eigen::Vector3d::Random();
    while (randVec == direction) randVec.setRandom();
    Eigen::Vector3d perp = direction.cross(randVec).normalized();
    newDirection = (direction + l*perp).normalized();
    
    return newDirection.dot(direction);
}

bool compare(const std::pair<double, double>& p1, const std::pair<double, double>& p2)
{
    return p1.first < p2.first;
}

double weightedAverage(std::vector<std::pair<double, double>>& coneDistancesAndWeights)
{
    int l = (int)coneDistancesAndWeights.size();
    if (l > 0) {
        // compute median
        std::sort(coneDistancesAndWeights.begin(), coneDistancesAndWeights.end(), compare);
        double median = 0;
        if (l % 2 == 0) median = coneDistancesAndWeights[l/2].first;
        else median = (coneDistancesAndWeights[l/2].first + coneDistancesAndWeights[l/2+1].first) / 2.0;
        
        // compute standard deviation
        double sd = 0;
        for (int i = 0; i < l; i++) {
            sd += pow(coneDistancesAndWeights[i].first-median, 2);
        }
        sd = sqrt(sd/(double)l);
        
        // compute weighted average
        double sum = 0;
        int count = 0;
        for (int i = 0; i < l; i++) {
            if (std::abs(coneDistancesAndWeights[i].first-median) < sd) {
                sum += coneDistancesAndWeights[i].first / coneDistancesAndWeights[i].second;
                count ++;
            }
        }
        
        return sum / (double)count;
    }
    
    return 0.0;
}

void Mesh::computeNormalizedSDFValues()
{
    Bvh bvh(this);
    double minSDF = INFINITY;
    double maxSDF = -INFINITY;
    Eigen::Vector3d direction;
    const int iterations = 30;
    
    // compute sdf as weighted average of cone rays for each face
    for (FaceIter f = faces.begin(); f != faces.end(); f++) {
        
        Eigen::Vector3d centroid = f->centroid();
        Eigen::Vector3d normal = -f->normal().normalized();
        
        std::vector<std::pair<double, double>> coneDistancesAndWeights;
        for (int i = 0; i < iterations; i++) {
            double coneAngle = randomDirectionInCone(direction, normal);
            double d = bvh.distance(centroid, direction);
            if (d != INFINITY) coneDistancesAndWeights.push_back(std::make_pair(d, coneAngle));
        }
        
        f->sdf = weightedAverage(coneDistancesAndWeights);
        if (f->sdf < minSDF) minSDF = f->sdf;
        if (f->sdf > maxSDF) maxSDF = f->sdf;
    }

    // normalize
    double log5 = log(5);
    double dSDF = maxSDF - minSDF;
    for (FaceIter f = faces.begin(); f != faces.end(); f++) {
        f->sdf = log(((f->sdf - minSDF) / dSDF)*4 + 1) / log5;
    }
}

void Mesh::cluster(const int k)
{
    // choose random centroids
    std::vector<double> centroids;
    for (int i = 0; i < k; i++) {
        int index = rand() % faces.size();
        while (std::find(centroids.begin(), centroids.end(), faces[index].sdf) != centroids.end()) {
            index = rand() % faces.size();
        }
        centroids.push_back(faces[index].sdf);
    }
    
    while (true) {
        // assign closest cluster to face
        for (FaceIter f = faces.begin(); f != faces.end(); f++) {
            double min = INFINITY;
            for (int i = 0; i < k; i++) {
                double distance = std::abs(f->sdf - centroids[i]);
                if (distance < min) {
                    min = distance;
                    f->cluster = i;
                }
            }
        }
        
        // update centroids
        std::vector<double> newCentroids;
        for (int i = 0; i < k; i++) {
            double sum = 0.0;
            int count = 0;
            for (FaceCIter f = faces.begin(); f != faces.end(); f++) {
                if (f->cluster == i) {
                    sum += f->sdf;
                    count ++;
                }
            }
            newCentroids.push_back(sum / (double)count);
        }
        
        // check for termination
        double distance = 0.0;
        for (int i = 0; i < k; i++) {
            distance += std::abs(newCentroids[i] - centroids[i]);
        }
        
        if (distance == 0.0) break;
        centroids = newCentroids;
    }
}

void Mesh::segment(const int k)
{
    // soft clustering
    cluster(k);
    
    // k graph cut
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
    }
    
    // determine radius
    double rMax = 0;
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
