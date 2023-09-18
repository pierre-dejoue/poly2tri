#pragma once

#include <poly2tri/common/shapes.h>
#include <poly2tri/sweep/cdt.h>

#include <vector>

bool TriangulationSanityChecks(const p2t::CDT& cdt, const std::vector<p2t::Triangle*>& triangles);

bool IsConstrainedDelaunay(const std::vector<p2t::Triangle*>& triangles);

std::vector<p2t::Point*> MakePointerVector(std::vector<p2t::Point>& points);
