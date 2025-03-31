#ifndef STRAINTENSOR_H_
#define STRAINTENSOR_H_

#include <vector>
#include <cmath>
#include "SystemStructures.h"
#include "System.h"

/*
 * Function: applyStrainToEdges
 * ----------------------------
 * Applies a strain tensor to update the rest lengths of edges.
 *
 * Parameters:
 *   epsilon_r       - Radial strain (e.g. +0.1 for 10% extension in radial direction)
 *   epsilon_t       - Tangential strain (can differ for anisotropic strain; same as epsilon_r for isotropic)
 *   nodeLocX        - x-coordinates of nodes
 *   nodeLocY        - y-coordinates of nodes
 *   edges2Nodes_1   - vector of indices for the first endpoint of each edge
 *   edges2Nodes_2   - vector of indices for the second endpoint of each edge
 *   centerX         - x-coordinate of the disc/tissue center
 *   centerY         - y-coordinate of the disc/tissue center
 *   edgeRestLengths - vector of edge rest lengths to be updated (in-place modification)
 */
void applyStrainToEdges(GeneralParams& generalParams, CoordInfoVecs& coordInfoVecs, LinearSpringInfoVecs& linearSpringInfoVecs);

//double epsilon_r, double epsilon_t,
                        //const std::vector<double>& nodeLocX,
                        //const std::vector<double>& nodeLocY,
                        //const std::vector<int>& edges2Nodes_1,
                        //const std::vector<int>& edges2Nodes_2,
                        //double centerX, double centerY,
                        //std::vector<double>& edgeRestLengths);

#endif // STRAINTENSOR_H_
