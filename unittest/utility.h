#pragma once

#include <vector>

namespace p2t
{
  class Triangle;
}

bool TriangulationSanityChecks(const std::vector<p2t::Triangle*>& triangles);

bool IsConstrainedDelaunay(const std::vector<p2t::Triangle*>& triangles);
