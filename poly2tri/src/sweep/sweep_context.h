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

#pragma once

#include <cassert>
#include <cstddef>
#include <vector>


namespace p2t {

// Inital triangle factor, seed triangle will extend 30% of
// PointSet width to both left and right.
constexpr double kAlpha = 0.3;

struct Point;
class Triangle;
struct Node;
struct Edge;

class SweepContext {
public:

  SweepContext();
  ~SweepContext();

  void AddPolyline(const Point* const* polyline, std::size_t num_points);
  void AddPolyline(const Point* polyline, std::size_t num_points, std::size_t stride = 0u);

  void AddHole(const Point* const* polyline, std::size_t num_points);
  void AddHole(const Point* polyline, std::size_t num_points, std::size_t stride = 0);

  void AddPoint(const Point* point);

  void AddPoints(const Point* const* points, std::size_t num_points);
  void AddPoints(const Point* points, std::size_t num_points, std::size_t stride = 0);

  const Point* head() const;

  const Point* tail() const;

  size_t point_count() const;

  void AddToMap(Triangle* triangle);

  const Point* GetPoint(size_t index);

  const std::vector<Edge*>& GetUpperEdges(size_t index);

  void MeshCleanExteriorTriangles(Triangle& interior_triangle);

  void MeshCleanHeadAndTail();

  const std::vector<Triangle*>& GetTriangles();

  const std::vector<Triangle*>& GetMap();

  void InitTriangulation();

  struct Basin {
    Node* left_node;
    Node* bottom_node;
    Node* right_node;
    double width;
    bool left_highest;

    Basin() :
      left_node(nullptr),
      bottom_node(nullptr),
      right_node(nullptr),
      width(0.0),
      left_highest(false)
    {
    }

    void Clear()
    {
      left_node = nullptr;
      bottom_node = nullptr;
      right_node = nullptr;
      width = 0.0;
      left_highest = false;
    }
  };

  struct EdgeEvent {
    Edge* constrained_edge;
    bool right;

    EdgeEvent() : constrained_edge(nullptr), right(false)
    { }
  };

private:

  template <typename GenPointPtr>
  void AddClosedPolylineGen(GenPointPtr generator, std::size_t num_points);

  template <typename GenPointPtr>
  void AddPointsGen(GenPointPtr generator, std::size_t num_points);

  friend class Sweep;

  struct SweepPoint
  {
    SweepPoint(const Point* p_) : p(p_), edges() { assert(p != nullptr); }

    const Point* p;
    std::vector<Edge*> edges;   // List of edges for which this point is the upper endpoint
  };
  static bool cmp(const SweepPoint& a, const SweepPoint& b);

  std::vector<SweepPoint> points_;
  std::vector<Edge*> edge_list_;
  std::vector<Triangle*> triangles_;
  std::vector<Triangle*> map_;

  // Artificial points added to the triangulation
  const Point* head_;
  const Point* tail_;

  Basin basin_;

  EdgeEvent edge_event_;

  void InitEdges(std::size_t polyline_begin_index, std::size_t num_points);

};

inline size_t SweepContext::point_count() const
{
  return points_.size();
}

inline const Point* SweepContext::head() const
{
  return head_;
}

inline const Point* SweepContext::tail() const
{
  return tail_;
}

}
