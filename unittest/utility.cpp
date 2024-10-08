#include "utility.h"

#include <cassert>
#include <iterator>
#include <numeric>
#include <set>

bool TriangulationSanityChecks(const p2t::CDT& cdt, const std::vector<p2t::Triangle*>& triangles)
{
  const auto all_input_points = cdt.GetInputPoints();
  std::set<const p2t::Point*> pointset(std::cbegin(all_input_points), std::cend(all_input_points));
  if (all_input_points.size() != pointset.size())
    return false;

  const auto& info = cdt.LastTriangulationInfo();
  if (static_cast<std::size_t>(info.nb_input_points) != pointset.size())
    return false;
  if (static_cast<std::size_t>(info.nb_output_triangles) != triangles.size())
    return false;

  std::set<p2t::Triangle*> triangleset(std::cbegin(triangles), std::cend(triangles));
  if (triangles.size() != triangleset.size())
    return false;

  for (p2t::Triangle* t : triangles) {
    if (t == nullptr)
      return false;
    for (int i = 0; i < 3; i++)
      if (pointset.count(t->GetPoint(i)) == 0)
        return false;
    for (int i = 0; i < 3; i++)
      if (t->GetNeighbor(i) && triangleset.count(t->GetNeighbor(i)) == 0)
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
      if (t->IsConstrainedEdge(e))
        continue;
      p2t::Triangle* ot = t->GetNeighbor(e);
      if (ot && (ot->GetState() == t->GetState())) {
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

std::size_t CountPoints(const std::vector<p2t::Triangle*>& triangles)
{
  std::set<const p2t::Point*> pointset;
  for (const auto* t: triangles)
  {
    pointset.insert(t->GetPoint(0));
    pointset.insert(t->GetPoint(1));
    pointset.insert(t->GetPoint(2));
  }
  return pointset.size();
}