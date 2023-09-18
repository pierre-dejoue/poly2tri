/*
 * Poly2Tri Copyright (c) 2009-2023, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <poly2tri/sweep/cdt.h>

#include "sweep_context.h"
#include "sweep.h"

#include <algorithm>
#include <iterator>

namespace p2t {

CDT::CDT() :
  sweep_context_(std::make_unique<SweepContext>())
{
}

CDT::CDT(const std::vector<Point*>& polyline) :
  CDT()
{
  AddPolyline(polyline);
}

CDT::~CDT() = default;

void CDT::AddPolyline(const std::vector<Point*>& polyline)
{
  sweep_context_->AddPolyline(polyline.data(), polyline.size());
}

void CDT::AddPolyline(const Point* const* polyline, std::size_t num_points)
{
  sweep_context_->AddPolyline(polyline, num_points);
}

void CDT::AddPolyline(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  sweep_context_->AddPolyline(polyline, num_points, stride);
}

void CDT::AddHole(const std::vector<Point*>& polyline)
{
  sweep_context_->AddHole(polyline.data(), polyline.size());
}

void CDT::AddHole(const Point* const* polyline, std::size_t num_points)
{
  sweep_context_->AddHole(polyline, num_points);
}

void CDT::AddHole(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  sweep_context_->AddHole(polyline, num_points, stride);
}

void CDT::AddPoint(const Point* point)
{
  sweep_context_->AddPoint(point);
}

void CDT::AddPoints(const std::vector<Point*>& points)
{
  sweep_context_->AddPoints(points.data(), points.size());
}

void CDT::AddPoints(const Point* const* points, std::size_t num_points)
{
  sweep_context_->AddPoints(points, num_points);
}

void CDT::AddPoints(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  sweep_context_->AddPoints(polyline, num_points, stride);
}

void CDT::Triangulate(Policy triangulation_policy)
{
  SweepContext& tcx = *sweep_context_;
  if (tcx.point_count() > 0) {
    tcx.InitTriangulation();
    Sweep(tcx).Triangulate(triangulation_policy);
  }
}

std::size_t CDT::GetTrianglesCount() const
{
  return sweep_context_->GetTriangles().size();
}

const std::vector<std::unique_ptr<Triangle>>& CDT::GetTriangles() const
{
  return sweep_context_->GetTriangles();
}

std::size_t CDT::GetInputPointsCount() const
{
  return sweep_context_->GetPoints().size();
}

std::vector<const Point*> CDT::GetInputPoints() const
{
  std::vector<const Point*> result;
  const auto& sweep_points = sweep_context_->GetPoints();
  result.reserve(sweep_points.size());
  std::transform(std::cbegin(sweep_points), std::cend(sweep_points), std::back_inserter(result), [](const auto& sweep_p) { return sweep_p.p; });
  return result;
}

std::size_t CDT::GetInputEdgesCount() const
{
  return sweep_context_->GetEdges().size();
}

std::vector<Edge*> CDT::GetInputEdges() const
{
  std::vector<Edge*> result;
  const auto& edges = sweep_context_->GetEdges();
  result.reserve(edges.size());
  std::transform(std::cbegin(edges), std::cend(edges), std::back_inserter(result), [](const auto& uniq_ptr) { return uniq_ptr.get(); });
  return result;
}

std::vector<Triangle*> GetTrianglesAsVector(const CDT& cdt)
{
  std::vector<p2t::Triangle*> result;
  result.reserve(cdt.GetTrianglesCount());
  cdt.GetTriangles(std::back_inserter(result));
  return result;
}

} // namespace p2t
