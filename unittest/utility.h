#pragma once

#include <vector>

namespace p2t
{
  class Triangle;
}

bool IsConstrainedDelaunay(const std::vector<p2t::Triangle*>& triangles);
