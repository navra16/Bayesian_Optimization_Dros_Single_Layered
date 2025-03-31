#include "StrainTensor.h"
#include <algorithm>   // for std::max
#include <cmath>
#include "System.h"

// This function updates each edge's rest length based on a strain tensor
// by decomposing the edge vector into radial and tangential components.
// The radial component is scaled by (1 + epsilon_r) and the tangential component
// by (1 + epsilon_t), and then they are recombined to form the new rest length.
void applyStrainToEdges( 
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs)

//double epsilon_r, double epsilon_t,
//                        const std::vector<double>& nodeLocX,
//                        const std::vector<double>& nodeLocY,
//                        const std::vector<int>& edges2Nodes_1,
//                        const std::vector<int>& edges2Nodes_2,
//                        double centerX, double centerY,
//                        std::vector<double>& edgeRestLengths)
{
    int numEdges = coordInfoVecs.edges2Nodes_1.size();
    for (int e = 0; e < numEdges; ++e) {
        int i = coordInfoVecs.edges2Nodes_1[e];
        int j = coordInfoVecs.edges2Nodes_2[e];
        double x1 = coordInfoVecs.nodeLocX[i], y1 = coordInfoVecs.nodeLocY[i];
        double x2 = coordInfoVecs.nodeLocX[j], y2 = coordInfoVecs.nodeLocY[j];

        // Compute the midpoint of the edge.
        double midX = 0.5 * (x1 + x2);
        double midY = 0.5 * (y1 + y2);

        // Compute the radial vector from the center to the midpoint.
        double rx = midX - generalParams.centerX;
        double ry = midY - generalParams.centerY;
        double r_norm = sqrt(rx * rx + ry * ry);
        double urx = (r_norm > 1e-12) ? rx / r_norm : 0.0;
        double ury = (r_norm > 1e-12) ? ry / r_norm : 0.0;

        // Compute the edge vector.
        double dx = x2 - x1;
        double dy = y2 - y1;
        double L = sqrt(dx * dx + dy * dy);

        // Decompose the edge vector into radial and tangential components.
        double L_r = dx * urx + dy * ury;
        // Use std::max to avoid negative values under the square root.
        double L_t = sqrt(std::max(0.0, L * L - L_r * L_r));

        // Apply the strain: scale radial and tangential components.
        double new_L_r = L_r * (1.0 + generalParams.epsilon_r);
        double new_L_t = L_t * (1.0 + generalParams.epsilon_t);

        // Recombine the components to compute the new rest length.
        double newRestLength = sqrt(new_L_r * new_L_r + new_L_t * new_L_t);

        // Update the edge's rest length in the host vector.
        linearSpringInfoVecs.edge_initial_length[e] = newRestLength;
    }
}
