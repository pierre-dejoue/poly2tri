#include "utility.h"

#include <poly2tri/common/shapes.h>

#include <cassert>

bool TriangulationSanityChecks(const std::vector<p2t::Triangle*>& triangles)
{
  for (p2t::Triangle* t : triangles) {
    if (t == nullptr)
      return false;
    if (!t->IsInterior())
      return false;
  }
  return true;
}

// Check that all edge that are not constrained ones are "locally optimal"
// which means that the quadrilateral formed by the union of the
// two adjacent triangles does not need to be legalized
bool IsConstrainedDelaunay(const std::vector<p2t::Triangle*>& triangles)
{
  for (p2t::Triangle* t : triangles) {
    // Loop on the triangle edges
    for (int e = 0; e < 3; e++) {
      if (t->constrained_edge[e])
        continue;
      p2t::Triangle* ot = t->GetNeighbor(e);
      if (ot && ot->IsInterior()) {
        const p2t::Point* p = t->GetPoint(e);
        const p2t::Point* op = ot->OppositePoint(*t, p);
        if (t->CircumcircleContains(*op)) {
          return false;
        }
      }
    }
  }
  return true;
}
