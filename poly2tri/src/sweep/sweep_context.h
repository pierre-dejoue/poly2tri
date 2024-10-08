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

#pragma once

#include <poly2tri/common/point.h>
#include <poly2tri/common/shapes.h>

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace p2t {

// Inital triangle factor, seed triangle will extend 30% of
// PointSet width to both left and right.
constexpr double kAlpha = 0.3;

struct Node;
struct Edge;

class TriangleStorage;
class NodeStorage;


// Wrapper structure around the Points provided by the user
struct SweepPoint
{
  struct UpperEdges
  {
    std::uint32_t nb_edges{0};
    std::array<std::uint32_t, 2> edges;
  };

  SweepPoint(const Point* p_) : p(p_), upper_edges{} { assert(p != nullptr); }
  SweepPoint(const SweepPoint&) = delete;
  SweepPoint& operator=(const SweepPoint&) = delete;
  SweepPoint(SweepPoint&&) noexcept = default;
  SweepPoint& operator=(SweepPoint&&) noexcept = default;

  static inline bool cmp(const SweepPoint& a, const SweepPoint& b) { return p2t::cmp(a.p, b.p); }

  const Point* p;
  UpperEdges upper_edges;         // List of constrained edges for which this point is the upper endpoint (max: 2)
};

// The sweep context
class SweepContext {
public:

  SweepContext();
  ~SweepContext();

  // Not moveable, not copyable
  SweepContext(SweepContext&&) = delete;
  SweepContext& operator=(SweepContext&&) = delete;
  SweepContext(const SweepContext&) = delete;
  SweepContext& operator=(const SweepContext&) = delete;

  void AddPolyline(const Point* const* polyline, std::size_t num_points, bool closed);
  void AddPolyline(const Point* polyline, std::size_t num_points, bool closed, std::size_t stride = 0u);

  void AddPoint(const Point* point);

  void AddPoints(const Point* const* points, std::size_t num_points);
  void AddPoints(const Point* points, std::size_t num_points, std::size_t stride = 0);

  void ClearInput();

  inline const Point* head() const { return &head_; }

  inline const Point* tail() const { return &tail_; }

  Triangle* AddTriangle(const Point* a, const Point* b, const Point* c);

  void DiscardTriangle(Triangle& t);

  Node* AddEmptyNode();

  Node* AddNode(const Point* p, Triangle* t = nullptr);

  void DiscardNode(Node* node);

  const Point* GetPoint(size_t index) const;

  const SweepPoint::UpperEdges& GetUpperEdges(size_t index) const;

  const Edge& GetConstrainedEdge(size_t index) const;

  void PopulateTriangleMap(Triangle::State_t filter);

  const std::vector<SweepPoint>& GetPoints() const;

  std::size_t GetEdgesCount() const { return constrained_edges_.size(); }

  const std::vector<Triangle*>& GetTriangles() const;

  void InitTriangulation();

  std::size_t InputMemoryFootprint() const;

  std::size_t NodeMemoryFootprint() const;

  std::size_t TriangleStorageFootprint() const;

  std::size_t TriangleStorageNbOfCreatedTriangles() const;

private:

  friend class Sweep;

  template <typename GenPointPtr>
  void AddPolylineGen(GenPointPtr generator, std::size_t num_points, bool closed);

  template <typename GenPointPtr>
  void AddPointsGen(GenPointPtr generator, std::size_t num_points);

  void InitEdges(std::size_t polyline_begin_index, std::size_t num_points, bool closed);

  void AllocateTriangleBuffer();

  // Triangulation input
  std::vector<SweepPoint> points_;
  std::vector<Edge>       constrained_edges_;

  // Triangle storage
  std::unique_ptr<TriangleStorage> triangle_storage_;

  // Node storage
  std::unique_ptr<NodeStorage> node_storage_;

  // The map of all triangles that are part of the triangulation
  std::vector<Triangle*> map_;

  // Artificial points added to the triangulation
  Point head_;
  Point tail_;

};

}
