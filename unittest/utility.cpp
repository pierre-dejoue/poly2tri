#include "utility.h"

#include <cassert>
#include <iterator>
#include <numeric>

bool TriangulationSanityChecks(const std::vector<p2t::Triangle*>& triangles)
{
  for (p2t::Triangle* t : triangles) {
    if (t == nullptr)
      return false;
    if (!t->IsInterior())
      return false;
    for (int i = 0; i < 3; i++)
      if (t->GetNeighbor(i) && !t->GetNeighbor(i)->IsInterior())
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

std::vector<p2t::Point*> MakePointerVector(std::vector<p2t::Point>& points)
{
  std::vector<p2t::Point*> result;
  result.resize(points.size(), nullptr);
  std::iota(result.begin(), result.end(), points.data());
  return result;
}
