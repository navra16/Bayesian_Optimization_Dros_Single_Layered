#ifndef SEGMENTVOLUMECONSTRAINT_H_
#define SEGMENTVOLUMECONSTRAINT_H_

#include "SystemStructures.h"
#include <vector>

// Applies a volume constraint to maintain the volume of a specified segment.
// The segment is defined by the indices of triangles and nodes that belong to it.
// The function computes the current volume from the triangle data, compares it to the equilibrium volume,
// and then applies a corrective force to the nodes in the segment.
void ApplySegmentVolumeConstraint(GeneralParams& generalParams,
                                  CoordInfoVecs& coordInfoVecs,
                                  const std::vector<int>& segmentTriangleIndices,
                                  const std::vector<int>& segmentNodeIndices);

#endif // SEGMENTVOLUMECONSTRAINT_H_
