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

#include "sweep_context.h"

#include <poly2tri/common/shapes.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <exception>
#include <iterator>
#include <numeric>
#include <stdexcept>


namespace p2t {

class TriangleStorage {
public:
  TriangleStorage(std::size_t nb_points)
  {
    //
    // Geometry tells us that the triangulation of a planer point set P of size n,
    // where k denotes the number of points in P that lie on the boundary of the
    // convex hull of P, has 2n−2−k triangles.
    // Assuming k is at least 3 in our use case we can evaluate the maximum
    // number of triangles required for the purpose of computing the CDT.
    //
    assert(nb_points >= 3);
    max_triangles_ = 2 * nb_points - 5;
    triangles_.reserve(max_triangles_);
  }

  std::size_t MemoryFootprint() const
  {
    return triangles_.capacity() * sizeof(Triangle);
  }

  Triangle* NewTriangle(const Point* a, const Point* b, const Point* c)
  {
    Triangle* result = nullptr;
    if (!discarded_triangles_.empty()) {
      // Re-use a discarded triangle
      result = discarded_triangles_.back();
      discarded_triangles_.pop_back();
      *result = Triangle(a, b, c);
    } else {
      // Initialize a Triangle from the storage (it was pre-allocated)
      assert(triangles_.capacity() >= max_triangles_);
      if (triangles_.size() >= max_triangles_) {
        throw std::runtime_error("TriangleStorage: Out of memory capacity");
      }
      result = &triangles_.emplace_back(a, b, c);
    }
    assert(result != nullptr);
    return result;
  }

  void DeleteTriangle(Triangle* t)
  {
    discarded_triangles_.push_back(t);
  }

private:
  std::size_t max_triangles_;
  std::vector<Triangle> triangles_;
  std::vector<Triangle*> discarded_triangles_;
};

SweepContext::SweepContext() = default;

SweepContext::~SweepContext() = default;

template <typename GenPointPtr>
void SweepContext::AddPolylineGen(GenPointPtr generator, std::size_t num_points, bool closed)
{
  if (closed && num_points < 3) {
    throw std::invalid_argument("A closed polyline must have at least 3 vertices");
  }
  if (!closed && num_points < 2) {
    throw std::invalid_argument("An open polyline must have at least 2 vertices");
  }
  const std::size_t begin_index = points_.size();
  points_.reserve(begin_index + num_points);
  std::generate_n(std::back_inserter(points_), num_points, [&generator]() { return SweepPoint(generator()); });
  InitEdges(begin_index, num_points, closed);
}

template <typename GenPointPtr>
void SweepContext::AddPointsGen(GenPointPtr generator, std::size_t num_points)
{
  const std::size_t begin_index = points_.size();
  points_.reserve(begin_index + num_points);
  std::generate_n(std::back_inserter(points_), num_points, [&generator]() { return SweepPoint(generator()); });
}

void SweepContext::AddPolyline(const Point* const* polyline, std::size_t num_points, bool closed)
{
  AddPolylineGen([&polyline]() { return *polyline++; }, num_points, closed);
}

void SweepContext::AddPolyline(const Point* polyline, std::size_t num_points, bool closed, std::size_t stride)
{
  if (stride == 0u) { stride = sizeof(Point); }
  auto mem_ptr = reinterpret_cast<const char*>(polyline);
  AddPolylineGen([&mem_ptr, stride]() { auto p = reinterpret_cast<const Point*>(mem_ptr); mem_ptr += stride; return p; }, num_points, closed);
}

void SweepContext::AddPoint(const Point* point)
{
  points_.emplace_back(point);
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

const std::vector<SweepPoint>& SweepContext::GetPoints() const
{
  return points_;
}

std::size_t SweepContext::GetEdgesCount() const
{
  return std::accumulate(std::cbegin(points_), std::cend(points_), std::size_t{0}, [](std::size_t cumul, const SweepPoint& point) {
    return cumul + point.edges.size();
  });
}

const std::vector<Triangle*>& SweepContext::GetTriangles() const
{
  return map_;
}

bool SweepPoint::cmp(const SweepPoint& a, const SweepPoint& b)
{
  return p2t::cmp(a.p, b.p);
}

void SweepContext::InitTriangulation()
{
  // Clear any previous triangulation
  map_.clear();
  triangle_storage_.reset();

  assert(!points_.empty());
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

  // Artificial points
  head_ = std::make_unique<const Point>(xmin - dx, ymin - dy);
  tail_ = std::make_unique<const Point>(xmax + dx, ymin - dy);

  // Sort input points along y-axis
  std::sort(points_.begin(), points_.end(), &SweepPoint::cmp);

  // Triangle buffer
  AllocateTriangleBuffer();
}

// Return the memory footprint of the Triangle storage, in bytes
std::size_t SweepContext::TriangleStorageFootprint() const
{
  return triangle_storage_ ? triangle_storage_->MemoryFootprint() : 0u;
}

void SweepContext::InitEdges(std::size_t polyline_begin_index, std::size_t num_points, bool closed)
{
  if (num_points == 0) {
    return;
  }
  const std::size_t begin = polyline_begin_index;
  const std::size_t end   = polyline_begin_index + num_points;
  const std::size_t last  = end - (closed ? 0 : 1);
  assert(end <= points_.size());
  for (std::size_t i = begin; i < last; i++) {
    std::size_t j = i < (end - 1) ? i + 1 : begin;
    Edge edge = Edge(points_[i].p, points_[j].p);
    const std::size_t upper_endpoint = (edge.q == points_[i].p ? i : j);
    points_[upper_endpoint].edges.emplace_back(std::move(edge));
  }
}

const Point* SweepContext::GetPoint(size_t index)
{
  return points_[index].p;
}

const std::vector<Edge>& SweepContext::GetUpperEdges(size_t index) const
{
  return points_[index].edges;
}

void SweepContext::AllocateTriangleBuffer()
{
  assert(points_.size() > 0);
  const std::size_t nb_points = points_.size() + 2;     // +2 artificial points: head and tail
  triangle_storage_ = std::make_unique<TriangleStorage>(nb_points);
}

Triangle* SweepContext::AddTriangleToMap(const Point* a, const Point* b, const Point* c)
{
  assert(triangle_storage_);
  Triangle* new_triangle = triangle_storage_->NewTriangle(a, b, c);
  map_.emplace_back(new_triangle);
  return new_triangle;
}

std::size_t SweepContext::MeshCleanExteriorTriangles()
{
  assert(triangle_storage_);
  const auto last_it = std::remove_if(std::begin(map_), std::end(map_), [this](auto& t) {
    if (!t->IsInterior()) {
      t->ClearNeighbors();    // To clear the reciprocate neighbor link
      triangle_storage_->DeleteTriangle(t);
      return true;
    }
    return false;
  });
  const auto erased = static_cast<std::size_t>(std::distance(last_it, std::end(map_)));
  map_.erase(last_it, std::end(map_));
  return erased;
}

} // namespace p2t
