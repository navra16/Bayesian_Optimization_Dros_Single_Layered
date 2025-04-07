#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"

// Declare the function ComputeLinearSprings, which calculates linear spring forces.
void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs);

// Define a functor for calculating linear spring forces for each edge.
struct LinearSpringFunctor {
    int SCALE_TYPE;
    bool nonuniform_wall_weakening_linear;
    double maxSpringScaler_linear;
    double scaling_pow;
    double gausssigma;
    double hilleqnconst;
    double hilleqnpow;
    double* scaling_per_edge;
    double spring_constant;
    double spring_constant_weak;
    // Remove the global rest length:
    // double length_zero;
    // New member: pointer to per-edge rest lengths:
    double* edge_rest_length;
    
    int* edges_in_upperhem;
    // int* boundaries_in_upperhem;

    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;

    // Constructor now accepts a pointer to the per-edge rest lengths.
    __host__ __device__ LinearSpringFunctor(
        int& _SCALE_TYPE,
        bool& _nonuniform_wall_weakening_linear,
        double& _maxSpringScaler_linear,
        double& _scaling_pow,
        double& _gausssigma,
        double& _hilleqnconst,
        double& _hilleqnpow,
        double* _scaling_per_edge,
        double& _spring_constant,
        double& _spring_constant_weak,
        /* Removed: double& _length_zero, */
        double* _edge_rest_length,           // New parameter
        int* _edges_in_upperhem,
        // int* _boundaries_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr)
      : SCALE_TYPE(_SCALE_TYPE),
        nonuniform_wall_weakening_linear(_nonuniform_wall_weakening_linear),
        maxSpringScaler_linear(_maxSpringScaler_linear),
        scaling_pow(_scaling_pow),
        gausssigma(_gausssigma),
        hilleqnconst(_hilleqnconst),
        hilleqnpow(_hilleqnpow),
        scaling_per_edge(_scaling_per_edge),
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        edge_rest_length(_edge_rest_length), // assign the pointer
        edges_in_upperhem(_edges_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

    // Functor operator: calculates force and energy for each spring.
    __device__ double operator()(const Tuuu &u3d) {
        int counter = thrust::get<0>(u3d);
        int place = 2 * counter;
        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

        if (edgeL != INT_MAX && edgeL >= 0 && edgeR != INT_MAX && edgeR >= 0) {
            double what_spring_constant;
            if (SCALE_TYPE == 0) {
                what_spring_constant = spring_constant * (1.0 - ((1.0/sqrt(2*3.14159*gausssigma)) *
                    exp(-(scaling_per_edge[counter]*scaling_per_edge[counter])/gausssigma)));
                if (what_spring_constant < spring_constant_weak) { what_spring_constant = spring_constant; }
            }
            else if (SCALE_TYPE == 1) {
                what_spring_constant = spring_constant_weak * pow(scaling_per_edge[counter], scaling_pow) +
                                       spring_constant_weak * (1 - pow(scaling_per_edge[counter], scaling_pow));
            }
            else if (SCALE_TYPE == 2) {
                what_spring_constant = spring_constant - (spring_constant - spring_constant_weak) * scaling_per_edge[counter];
            }
            else if (SCALE_TYPE == 3) {
                if (edges_in_upperhem[counter] == 1) {
                    what_spring_constant = spring_constant_weak;
                }
                else if (edges_in_upperhem[counter] == 0) {
                    what_spring_constant = (spring_constant_weak + spring_constant) / 2.0;
                }
                else {
                    what_spring_constant = spring_constant;
                }
            }
            else if (SCALE_TYPE == 4) {
                if (nonuniform_wall_weakening_linear == true) {
                    double spectrum = maxSpringScaler_linear * spring_constant - spring_constant_weak;
                    what_spring_constant = spring_constant_weak + ((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[counter], hilleqnpow))) * spectrum);
                    if (what_spring_constant < spring_constant_weak) { what_spring_constant = spring_constant_weak; }
                }
                else {
                    if (edges_in_upperhem[counter] == 1) {
                        what_spring_constant = spring_constant_weak;
                    }
                    else if (edges_in_upperhem[counter] == 0) {
                        what_spring_constant = (spring_constant_weak + spring_constant) / 2.0;
                    }
                    else {
                        what_spring_constant = spring_constant;
                    }
                }
            }

            // Compute forces between node edgeL and edgeR.
            double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
            double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
            double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];
            double length_current = sqrt(xLoc_LR*xLoc_LR + yLoc_LR*yLoc_LR + zLoc_LR*zLoc_LR);
            double energy = 0.0;
            //edge_rest_length[counter] = length_current;
            // Here, instead of using the single global rest length "length_zero", we use the per-spring value:
            double restLen = edge_rest_length[counter];

            if (length_current != restLen) {
                double magnitude = -(what_spring_constant) * (length_current - restLen);
                
                idKey[place] = edgeL;
                forceXAddr[place] = magnitude * (xLoc_LR / length_current);
                forceYAddr[place] = magnitude * (yLoc_LR / length_current);
                forceZAddr[place] = magnitude * (zLoc_LR / length_current);

                idKey[place + 1] = edgeR;
                forceXAddr[place + 1] = -magnitude * (xLoc_LR / length_current);
                forceYAddr[place + 1] = -magnitude * (yLoc_LR / length_current);
                forceZAddr[place + 1] = -magnitude * (zLoc_LR / length_current);
                
                energy = (what_spring_constant / 2.0) * (length_current - restLen) * (length_current - restLen);
            }
            return energy;
        } else {
            return 0.0;
        }
    }
};

#endif
