#include "SegmentVolumeConstraint.h"
#include <cmath>
#include <iostream>

// Helper function: Compute the volume of a closed surface using the divergence theorem,
// summing contributions from each triangle in the segment.
// Here we assume each triangle contributes: V_triangle = (1/6) * (x1*(y2*z3 - y3*z2) + ...).
// For a non-closed surface, this is an approximation.
double computeSegmentVolume(const CoordInfoVecs& coordInfoVecs,
                            const std::vector<int>& segmentTriangleIndices) {
    double volume = 0.0;
    int numTriangles = segmentTriangleIndices.size();
    for (int idx = 0; idx < numTriangles; ++idx) {
        int t = segmentTriangleIndices[idx];
        int i = coordInfoVecs.triangles2Nodes_1[t];
        int j = coordInfoVecs.triangles2Nodes_2[t];
        int k = coordInfoVecs.triangles2Nodes_3[t];
        double x1 = coordInfoVecs.nodeLocX[i], y1 = coordInfoVecs.nodeLocY[i], z1 = coordInfoVecs.nodeLocZ[i];
        double x2 = coordInfoVecs.nodeLocX[j], y2 = coordInfoVecs.nodeLocY[j], z2 = coordInfoVecs.nodeLocZ[j];
        double x3 = coordInfoVecs.nodeLocX[k], y3 = coordInfoVecs.nodeLocY[k], z3 = coordInfoVecs.nodeLocZ[k];
        double vol = (x1 * (y2 * z3 - y3 * z2) +
                      x2 * (y3 * z1 - y1 * z3) +
                      x3 * (y1 * z2 - y2 * z1)) / 6.0;
        volume += vol;
    }
    return fabs(volume);
}

void ApplySegmentVolumeConstraint(GeneralParams& generalParams,
                                  CoordInfoVecs& coordInfoVecs,
                                  const std::vector<int>& segmentTriangleIndices,
                                  const std::vector<int>& segmentNodeIndices)
{
    // Compute the current volume of the segment.
    double currentVolume = computeSegmentVolume(coordInfoVecs, segmentTriangleIndices);
    // Use the equilibrium volume from generalParams (or define a new parameter if desired).
    double eqVolume = generalParams.eq_total_volume;  
    // Compute volume error.
    double volumeError = currentVolume - eqVolume;
    
    // Use a volume spring constant (assumed stored in generalParams.volume_spring_constant) to determine the corrective force.
    double K_volume_segment = generalParams.volume_spring_constant;
    double correctiveForce = K_volume_segment * volumeError;
    
    // For simplicity, apply the corrective force uniformly to all nodes in the segment in the Z direction.
    int numNodes = segmentNodeIndices.size();
    if (numNodes == 0) return;
    double forcePerNode = correctiveForce / numNodes;
    
    // Update the Z-component of the node force for each node in the segment.
    for (int n : segmentNodeIndices) {
        // Note: this assumes that coordInfoVecs.nodeForceZ is accessible on host.
        coordInfoVecs.nodeForceZ[n] -= forcePerNode;
    }
    
    std::cout << "Segment volume: " << currentVolume 
              << " (eq: " << eqVolume << "), applied corrective force (Z): " << correctiveForce << std::endl;
}
