/*
 * Poly2Tri Copyright (c) 2009-2023, Poly2Tri Contributors
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
#include <exception>
#include <stdexcept>
#include <iterator>

namespace p2t {

namespace {
  constexpr bool CLOSED_POLYLINE = true;
  constexpr bool OPEN_POLYLINE = false;
}

CDT::CDT() :
  sweep_context_(std::make_unique<SweepContext>()),
  info_{}
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
  if (!sweep_context_->GetPoints().empty()) {
    throw std::invalid_argument("The outer polyline must be added first and only once");
  }
  sweep_context_->AddPolyline(polyline.data(), polyline.size(), CLOSED_POLYLINE);
}

void CDT::AddPolyline(const Point* const* polyline, std::size_t num_points)
{
  if (!sweep_context_->GetPoints().empty()) {
    throw std::invalid_argument("The outer polyline must be added first and only once");
  }
  sweep_context_->AddPolyline(polyline, num_points, CLOSED_POLYLINE);
}

void CDT::AddPolyline(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  if (!sweep_context_->GetPoints().empty()) {
    throw std::invalid_argument("The outer polyline must be added first and only once");
  }
  sweep_context_->AddPolyline(polyline, num_points, CLOSED_POLYLINE, stride);
}

void CDT::AddHole(const std::vector<Point*>& polyline)
{
  sweep_context_->AddPolyline(polyline.data(), polyline.size(), CLOSED_POLYLINE);
}

void CDT::AddHole(const Point* const* polyline, std::size_t num_points)
{
  sweep_context_->AddPolyline(polyline, num_points, CLOSED_POLYLINE);
}

void CDT::AddHole(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  sweep_context_->AddPolyline(polyline, num_points, CLOSED_POLYLINE, stride);
}

void CDT::AddOpenPolyline(const Point* const* polyline, std::size_t num_points)
{
  sweep_context_->AddPolyline(polyline, num_points, OPEN_POLYLINE);
}

void CDT::AddOpenPolyline(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  sweep_context_->AddPolyline(polyline, num_points, OPEN_POLYLINE, stride);
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
  info_ = Info{};
  if (!tcx.GetPoints().empty()) {
    tcx.InitTriangulation();
    info_.nb_input_points = static_cast<unsigned int>(tcx.GetPoints().size());
    info_.nb_input_edges = static_cast<unsigned int>(tcx.GetEdgesCount());
    info_.triangles_memory_footprint_in_bytes = tcx.TriangleStorageFootprint();
    Sweep(tcx, info_).Triangulate(triangulation_policy);
  }
}

std::size_t CDT::GetTrianglesCount() const
{
  return sweep_context_->GetTriangles().size();
}

const std::vector<Triangle*>& CDT::GetTriangles() const
{
  return sweep_context_->GetTriangles();
}

const CDT::Info& CDT::LastTriangulationInfo() const
{
  return info_;
}

std::vector<const Point*> CDT::GetInputPoints() const
{
  std::vector<const Point*> result;
  const auto& sweep_points = sweep_context_->GetPoints();
  result.reserve(sweep_points.size());
  std::transform(std::cbegin(sweep_points), std::cend(sweep_points), std::back_inserter(result), [](const auto& sweep_p) { return sweep_p.p; });
  return result;
}

} // namespace p2t
