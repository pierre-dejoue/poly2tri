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

#include "sweep_context.h"

#include <poly2tri/common/shapes.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <exception>
#include <iterator>
#include <stdexcept>


namespace p2t {

SweepContext::SweepContext() :
  points_(),
  edge_list_(),
  triangles_(),
  map_(),
  head_(nullptr),
  tail_(nullptr),
  basin_(),
  edge_event_()
{
}

template <typename GenPointPtr>
void SweepContext::AddClosedPolylineGen(GenPointPtr generator, std::size_t num_points)
{
  if (num_points < 3) {
    throw std::invalid_argument("A polyline must have at least 3 vertices");
  }
  const std::size_t begin_index = points_.size();
  points_.reserve(begin_index + num_points);
  std::generate_n(std::back_inserter(points_), num_points, [&generator]() { return SweepPoint(generator()); });
  InitEdges(begin_index, num_points);
}

void SweepContext::AddPolyline(const Point* const* polyline, std::size_t num_points)
{
  if (!points_.empty()) {
    throw std::invalid_argument("The outer polyline must be added first and only once");
  }
  AddClosedPolylineGen([&polyline]() { return *polyline++; }, num_points);
}

void SweepContext::AddPolyline(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  if (!points_.empty()) {
    throw std::invalid_argument("The outer polyline must be added first and only once");
  }
  if (stride == 0u) { stride = sizeof(Point); }
  auto mem_ptr = reinterpret_cast<const char*>(polyline);
  AddClosedPolylineGen([&mem_ptr, stride]() { auto p = reinterpret_cast<const Point*>(mem_ptr); mem_ptr += stride; return p; }, num_points);
}

void SweepContext::AddHole(const Point* const* polyline, std::size_t num_points)
{
  AddClosedPolylineGen([&polyline]() { return *polyline++; }, num_points);
}

void SweepContext::AddHole(const Point* polyline, std::size_t num_points, std::size_t stride)
{
  if (stride == 0u) { stride = sizeof(Point); }
  auto mem_ptr = reinterpret_cast<const char*>(polyline);
  AddClosedPolylineGen([&mem_ptr, stride]() { auto p = reinterpret_cast<const Point*>(mem_ptr); mem_ptr += stride; return p; }, num_points);
}

void SweepContext::AddPoint(const Point* point)
{
  points_.emplace_back(point);
}

template <typename GenPointPtr>
void SweepContext::AddPointsGen(GenPointPtr generator, std::size_t num_points)
{
  const std::size_t begin_index = points_.size();
  points_.reserve(begin_index + num_points);
  std::generate_n(std::back_inserter(points_), num_points, [&generator]() { return SweepPoint(generator()); });
}

void SweepContext::AddPoints(const Point* const* points, std::size_t num_points)
{
  AddPointsGen([&points]() { return *points++; }, num_points);
}

void SweepContext::AddPoints(const Point* points, std::size_t num_points, std::size_t stride)
{
  if (stride == 0u) { stride = sizeof(Point); }
  auto mem_ptr = reinterpret_cast<const char*>(points);
  AddPointsGen([&mem_ptr, stride]() { auto p = reinterpret_cast<const Point*>(mem_ptr); mem_ptr += stride; return p; }, num_points);
}

const std::vector<Triangle*>& SweepContext::GetTriangles()
{
  return triangles_;
}

const std::vector<Triangle*>& SweepContext::GetMap()
{
  return map_;
}


bool SweepContext::cmp(const SweepPoint& a, const SweepPoint& b)
{
  return p2t::cmp(a.p, b.p);
}

void SweepContext::InitTriangulation()
{
  assert(points_.size() > 0);
  double xmax(points_[0].p->x), xmin(points_[0].p->x);
  double ymax(points_[0].p->y), ymin(points_[0].p->y);

  // Calculate bounds
  for (const auto& sweep_point : points_) {
    const Point& p = *sweep_point.p;
    if (p.x > xmax)
      xmax = p.x;
    if (p.x < xmin)
      xmin = p.x;
    if (p.y > ymax)
      ymax = p.y;
    if (p.y < ymin)
      ymin = p.y;
  }

  double dx = kAlpha * (xmax - xmin);
  double dy = kAlpha * (ymax - ymin);
  head_ = new Point(xmin - dx, ymin - dy);
  tail_ = new Point(xmax + dx, ymin - dy);

  // Sort points along y-axis
  std::sort(points_.begin(), points_.end(), &SweepContext::cmp);
}

void SweepContext::InitEdges(std::size_t polyline_begin_index, std::size_t num_points)
{
  assert(num_points > 1);
  const std::size_t begin = polyline_begin_index;
  const std::size_t end   = polyline_begin_index + num_points;
  assert(end <= points_.size());
  for (std::size_t i = begin; i < end; i++) {
    std::size_t j = i < (end - 1) ? i + 1 : begin;
    edge_list_.emplace_back(new Edge(points_[i].p, points_[j].p));
    std::size_t upper_endpoint = (edge_list_.back()->q == points_[i].p ? i : j);
    points_[upper_endpoint].edges.emplace_back(edge_list_.back());
  }
}

const Point* SweepContext::GetPoint(size_t index)
{
  return points_[index].p;
}

const std::vector<Edge*>& SweepContext::GetUpperEdges(size_t index)
{
  return points_[index].edges;
}


void SweepContext::AddToMap(Triangle* triangle)
{
  map_.push_back(triangle);
}

void SweepContext::MeshCleanExteriorTriangles(Triangle& interior_triangle)
{
  std::vector<Triangle *> triangles;
  triangles.push_back(&interior_triangle);

  while(!triangles.empty()) {
    Triangle *t = triangles.back();
    triangles.pop_back();

    if (t != nullptr && !t->IsInterior()) {
      t->IsInterior(true);
      triangles_.push_back(t);
      for (int i = 0; i < 3; i++) {
        if (!t->constrained_edge[i])
          triangles.push_back(t->GetNeighbor(i));
      }
    }
  }
}

void SweepContext::MeshCleanHeadAndTail()
{
  assert(triangles_.empty());
  for (Triangle* t : map_)
  {
    const Point* p0 =  t->GetPoint(0);
    const Point* p1 =  t->GetPoint(1);
    const Point* p2 =  t->GetPoint(2);
    if (p0 != head_ && p0 != tail_ &&
        p1 != head_ && p1 != tail_ &&
        p2 != head_ && p2 != tail_)
    {
      t->IsInterior(true);
      triangles_.emplace_back(t);
    }
  }
}

SweepContext::~SweepContext()
{
    // Clean up memory
    delete head_;
    delete tail_;

    for (auto ptr : map_) {
      delete ptr;
    }
    // pointers in triangles_ are copies of pointers in map_ and don't need to be deleted

    for (auto& i : edge_list_) {
      delete i;
    }
}

} // namespace p2t
