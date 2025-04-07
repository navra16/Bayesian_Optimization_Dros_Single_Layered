#ifndef STRAINTENSOR_H_
#define STRAINTENSOR_H_

#include <thrust/device_vector.h>
#include "SystemStructures.h"
#include "System.h"

// Functor to update a single edge’s rest length based on a spatially varying strain tensor.
// The radial strain is computed by a linear interpolation between epsilon_r_center and epsilon_r_edge
// according to the distance of the edge midpoint from the disc center. The tangential strain (epsilon_t)
// is applied uniformly.
struct StrainTensorFunctor {
    double epsilon_r_center;   // radial strain at the center
    double epsilon_r_edge;     // radial strain at the boundary
    double epsilon_t;          // uniform tangential strain
    double centerX;            // disc center X
    double centerY;            // disc center Y
    double discRadius;         // disc radius

    // Pointers to device data for node positions and edge connectivity.
    const double* nodeLocX;
    const double* nodeLocY;
    const int* edges2Nodes1;
    const int* edges2Nodes2;
    // Pointer to the device vector of rest lengths to update.
    double* edgeRest;

    __host__ __device__
    StrainTensorFunctor(double _epsilon_r_center, double _epsilon_r_edge, double _epsilon_t,
                        double _centerX, double _centerY, double _discRadius,
                        const double* _nodeLocX, const double* _nodeLocY,
                        const int* _edges2Nodes1, const int* _edges2Nodes2,
                        double* _edgeRest)
      : epsilon_r_center(_epsilon_r_center), epsilon_r_edge(_epsilon_r_edge),
        epsilon_t(_epsilon_t), centerX(_centerX), centerY(_centerY), discRadius(_discRadius),
        nodeLocX(_nodeLocX), nodeLocY(_nodeLocY),
        edges2Nodes1(_edges2Nodes1), edges2Nodes2(_edges2Nodes2),
        edgeRest(_edgeRest) {}

    __device__
    void operator()(const int edgeID) const {
        int i = edges2Nodes1[edgeID];
        int j = edges2Nodes2[edgeID];

        double x1 = nodeLocX[i];
        double y1 = nodeLocY[i];
        double x2 = nodeLocX[j];
        double y2 = nodeLocY[j];

        // Compute the midpoint of the edge.
        double midX = 0.5 * (x1 + x2);
        double midY = 0.5 * (y1 + y2);

        // Compute the radial vector from the disc center to the midpoint.
        double rx = midX - centerX;
        double ry = midY - centerY;
        double r_norm = sqrt(rx * rx + ry * ry);
        double urx = (r_norm > 1e-12) ? rx / r_norm : 0.0;
        double ury = (r_norm > 1e-12) ? ry / r_norm : 0.0;

        // Compute the edge vector.
        double dx = x2 - x1;
        double dy = y2 - y1;
        double L = sqrt(dx * dx + dy * dy);

        // Decompose the edge vector into radial and tangential components.
        double L_r = dx * urx + dy * ury;
        double L_t = sqrt(fmax(0.0, L * L - L_r * L_r));

        // Compute local radial strain: linearly interpolate between center and edge strains.
        double fraction = (discRadius > 1e-12) ? (r_norm / discRadius) : 0.0;
        if (fraction > 1.0) fraction = 1.0;
        double local_epsilon_r = epsilon_r_center * (1.0 - fraction) + epsilon_r_edge * fraction;

        // Apply the strain: scale radial and tangential components.
        double new_L_r = L_r * (1.0 + local_epsilon_r);
        double new_L_t = L_t * (1.0 + epsilon_t);
        double newRestLength = sqrt(new_L_r * new_L_r + new_L_t * new_L_t);

        // Update the edge's rest length.
        edgeRest[edgeID] = newRestLength;
    }
};

// Function declaration that applies the strain tensor using the functor.
void applyStrainToEdges(GeneralParams& generalParams,
                        CoordInfoVecs& coordInfoVecs,
                        LinearSpringInfoVecs& linearSpringInfoVecs);

#endif // STRAINTENSOR_H_


//#ifndef STRAINTENSOR_H_
//#define STRAINTENSOR_H_
//
//#include <thrust/device_vector.h>
//#include "SystemStructures.h"
//#include "System.h"
//
//// Functor to update a single edge’s rest length based on the strain tensor.
//struct StrainTensorFunctor {
//    double epsilon_r;   // radial strain increment
//    double epsilon_t;   // tangential strain increment
//    double centerX;     // disc center X
//    double centerY;     // disc center Y
//
//    // Pointers to device data for node positions and edge connectivity.
//    const double* nodeLocX;
//    const double* nodeLocY;
//    const int* edges2Nodes1;
//    const int* edges2Nodes2;
//    // Pointer to the device vector of rest lengths to update.
//    double* edgeRest;
//
//    __host__ __device__
//    StrainTensorFunctor(double _epsilon_r, double _epsilon_t,
//                        double _centerX, double _centerY,
//                        const double* _nodeLocX, const double* _nodeLocY,
//                        const int* _edges2Nodes1, const int* _edges2Nodes2,
//                        double* _edgeRest)
//      : epsilon_r(_epsilon_r), epsilon_t(_epsilon_t),
//        centerX(_centerX), centerY(_centerY),
//        nodeLocX(_nodeLocX), nodeLocY(_nodeLocY),
//        edges2Nodes1(_edges2Nodes1), edges2Nodes2(_edges2Nodes2),
//        edgeRest(_edgeRest) {}
//
//    __device__
//    void operator()(const int edgeID) const {
//        int i = edges2Nodes1[edgeID];
//        int j = edges2Nodes2[edgeID];
//
//        double x1 = nodeLocX[i];
//        double y1 = nodeLocY[i];
//        double x2 = nodeLocX[j];
//        double y2 = nodeLocY[j];
//
//        double midX = 0.5 * (x1 + x2);
//        double midY = 0.5 * (y1 + y2);
//
//        double rx = midX - centerX;
//        double ry = midY - centerY;
//        double r_norm = sqrt(rx*rx + ry*ry);
//        double urx = (r_norm > 1e-12) ? rx/r_norm : 0.0;
//        double ury = (r_norm > 1e-12) ? ry/r_norm : 0.0;
//
//        double dx = x2 - x1;
//        double dy = y2 - y1;
//        double L = sqrt(dx*dx + dy*dy);
//
//        double L_r = dx * urx + dy * ury;
//        double L_t = sqrt(fmax(0.0, L*L - L_r*L_r));
//
//        double new_L_r = L_r * (1.0 + epsilon_r);
//        double new_L_t = L_t * (1.0 + epsilon_t);
//        double newRestLength = sqrt(new_L_r*new_L_r + new_L_t*new_L_t);
//
//        edgeRest[edgeID] = newRestLength;
//    }
//};
//
//void applyStrainToEdges(GeneralParams& generalParams,
//                        CoordInfoVecs& coordInfoVecs,
//                        LinearSpringInfoVecs& linearSpringInfoVecs);
//
//#endif // STRAINTENSOR_H_
//
//
////#ifndef STRAINTENSOR_H_
////#define STRAINTENSOR_H_
////
////#include <vector>
////#include <cmath>
////#include "SystemStructures.h"
////#include "System.h"
////
/////*
//// * Function: applyStrainToEdges
//// * ----------------------------
//// * Applies a strain tensor to update the rest lengths of edges.
//// *
//// * Parameters:
//// *   epsilon_r       - Radial strain (e.g. +0.1 for 10% extension in radial direction)
//// *   epsilon_t       - Tangential strain (can differ for anisotropic strain; same as epsilon_r for isotropic)
//// *   nodeLocX        - x-coordinates of nodes
//// *   nodeLocY        - y-coordinates of nodes
//// *   edges2Nodes_1   - vector of indices for the first endpoint of each edge
//// *   edges2Nodes_2   - vector of indices for the second endpoint of each edge
//// *   centerX         - x-coordinate of the disc/tissue center
//// *   centerY         - y-coordinate of the disc/tissue center
//// *   edgeRestLengths - vector of edge rest lengths to be updated (in-place modification)
//// */
////void applyStrainToEdges(GeneralParams& generalParams, CoordInfoVecs& coordInfoVecs, LinearSpringInfoVecs& linearSpringInfoVecs);
////
//////double epsilon_r, double epsilon_t,
////                        //const std::vector<double>& nodeLocX,
////                        //const std::vector<double>& nodeLocY,
////                        //const std::vector<int>& edges2Nodes_1,
////                        //const std::vector<int>& edges2Nodes_2,
////                        //double centerX, double centerY,
////                        //std::vector<double>& edgeRestLengths);
////
////#endif // STRAINTENSOR_H_
