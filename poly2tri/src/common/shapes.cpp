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
#include <poly2tri/common/shapes.h>
#include "utils.h"

#include <cassert>
#include <iostream>
#include <utility>

namespace p2t {

std::ostream& operator<<(std::ostream& out, const Point& point)
{
  return out << "{ " << point.x << ", " << point.y << " }";
}

Edge::Edge(const Point* p1, const Point* p2) : p(p1), q(p2)
{
  assert(p1); assert(p2);
  if (p1->y > p2->y) {
    std::swap(p, q);
  } else if (p1->y == p2->y) {
    if (p1->x > p2->x) {
      std::swap(p, q);
    } else if (p1->x == p2->x) {
      // Repeat points
      throw std::runtime_error("Edge::Edge: *p1 == *p2");
    }
  }
}

std::ostream& operator<<(std::ostream& out, const Edge& edge)
{
  return out << "{ p=" << *edge.p << ", q=" << *edge.q << " }";
}

std::ostream& operator<<(std::ostream& out, const Triangle& t)
{
  return out << "{ a=" << *t.GetPoint(0) << ", b=" << *t.GetPoint(1) << ", c=" << *t.GetPoint(2) << " }";
}

std::ostream& operator<<(std::ostream& out, Orientation o)
{
  switch (o)
  {
    case Orientation::COLLINEAR:
      return out << "COLLINEAR";
    case Orientation::CCW:
      return out << "CCW";
    case Orientation::CW:
      return out << "CW";
    default:
      assert(0);
      return out;
  }
}

Triangle::Triangle(const Point* a, const Point* b, const Point* c)
{
  assert(a); assert(b); assert(c);
  points_[0] = a; points_[1] = b; points_[2] = c;
  neighbors_[0] = nullptr; neighbors_[1] = nullptr; neighbors_[2] = nullptr;
  constrained_edge_[0] = constrained_edge_[1] = constrained_edge_[2] = false;
  delaunay_edge_[0] = delaunay_edge_[1] = delaunay_edge_[2] = false;
  interior_ = false;
  assert(Orient2d(*points_[0], *points_[1], *points_[2]) != Orientation::CW);
}

// Update neighbor pointers
void Triangle::MarkNeighbor(const Point* p1, const Point* p2, Triangle* t)
{
  if ((p1 == points_[2] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[2]))
    neighbors_[0] = t;
  else if ((p1 == points_[0] && p2 == points_[2]) || (p1 == points_[2] && p2 == points_[0]))
    neighbors_[1] = t;
  else if ((p1 == points_[0] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[0]))
    neighbors_[2] = t;
  else
    assert(0);
}

// Exhaustive search to update neighbor pointers
void Triangle::MarkNeighbor(Triangle& t)
{
  if (t.Contains(points_[1], points_[2])) {
    neighbors_[0] = &t;
    t.MarkNeighbor(points_[1], points_[2], this);
  } else if (t.Contains(points_[0], points_[2])) {
    neighbors_[1] = &t;
    t.MarkNeighbor(points_[0], points_[2], this);
  } else if (t.Contains(points_[0], points_[1])) {
    neighbors_[2] = &t;
    t.MarkNeighbor(points_[0], points_[1], this);
  } else {
    assert(0);
  }
}

// The caller MUST clear the link in the other direction
void Triangle::ClearNeighbor(const Triangle* triangle)
{
  if (neighbors_[0] == triangle) {
    neighbors_[0] = nullptr;
    delaunay_edge_[0] = false;
  } else if (neighbors_[1] == triangle) {
    neighbors_[1] = nullptr;
    delaunay_edge_[1] = false;
  } else {
    assert(neighbors_[2] == triangle);
    neighbors_[2] = nullptr;
    delaunay_edge_[2] = false;
  }
}

// Disconnect this triangle from all its neighbors
void Triangle::ClearNeighbors()
{
  if (neighbors_[0]) { neighbors_[0]->ClearNeighbor(this); }
  neighbors_[0] = nullptr;
  delaunay_edge_[0] = false;
  if (neighbors_[1]) { neighbors_[1]->ClearNeighbor(this); }
  neighbors_[1] = nullptr;
  delaunay_edge_[1] = false;
  if (neighbors_[2]) { neighbors_[2]->ClearNeighbor(this); }
  neighbors_[2] = nullptr;
  delaunay_edge_[2] = false;
}

void Triangle::ClearDelaunayEdges()
{
  delaunay_edge_[0] = delaunay_edge_[1] = delaunay_edge_[2] = false;
}

const Point* Triangle::OppositePoint(Triangle& t, const Point* p)
{
  const Point* cw = t.PointCW(p);
  return PointCW(cw);
}

// Legalize triangle by rotating clockwise around point
// The indices of the opposing points after the legalization are the same as before (although those are not the same points)
void Triangle::Legalize(const Point* p, const Point* op)
{
  assert(!IsConstrainedEdge(p));

  if (p == points_[0]) {
    points_[0] = points_[2];
    // neighbors_[0] remains the same (it is the opposite triangle 'ot')
    // constrained_edge_[0] is false;
    // delaunay_edge_[0] remains the same

    points_[2] = op;
    neighbors_[2] = neighbors_[1];
    constrained_edge_[2] = constrained_edge_[1];
    // delaunay_edge_[2] set below after all points are shifted

    points_[1] = p;
    neighbors_[1] = nullptr;      // To be set by the caller
    // constrained_edge[1] to be set by the caller
    // delaunay_edge[1] to be set false by the caller

    SetDelaunayEdge(2, false);

  } else if (p == points_[1]) {
    points_[1] = points_[0];
    // neighbors_[1] remains the same (it is the opposite triangle 'ot')
    // constrained_edge_[1] is false;
    // delaunay_edge_[1] remains the same

    points_[0] = op;
    neighbors_[0] = neighbors_[2];
    constrained_edge_[0] = constrained_edge_[2];
    // delaunay_edge_[0] set below after all points are shifted

    points_[2] = p;
    neighbors_[2] = nullptr;      // To be set by the caller
    // constrained_edge[2] to be set by the caller
    // delaunay_edge[2] to be set false by the caller

    SetDelaunayEdge(0, false);

  } else {
    assert(p == points_[2]);
    points_[2] = points_[1];
    // neighbors_[2] remains the same (it is the opposite triangle 'ot')
    // constrained_edge_[2] is false;
    // delaunay_edge_[2] remains the same

    points_[1] = op;
    neighbors_[1] = neighbors_[0];
    constrained_edge_[1] = constrained_edge_[0];
    // delaunay_edge_[1] set below after all points are shifted

    points_[0] = p;
    neighbors_[0] = nullptr;      // To be set by the caller
    // constrained_edge[0] to be set by the caller
    // delaunay_edge[0] to be set false by the caller

    SetDelaunayEdge(1, false);
  }
}

int Triangle::Index(const Point* p)
{
  if (p == points_[0]) {
    return 0;
  } else if (p == points_[1]) {
    return 1;
  } else if (p == points_[2]) {
    return 2;
  }
  assert(0);
  return -1;
}

int Triangle::EdgeIndex(const Point* p1, const Point* p2)
{
  if (points_[0] == p1) {
    if (points_[1] == p2) {
      return 2;
    } else if (points_[2] == p2) {
      return 1;
    }
  } else if (points_[1] == p1) {
    if (points_[2] == p2) {
      return 0;
    } else if (points_[0] == p2) {
      return 2;
    }
  } else if (points_[2] == p1) {
    if (points_[0] == p2) {
      return 1;
    } else if (points_[1] == p2) {
      return 0;
    }
  }
  return -1;
}

void Triangle::SetConstrainedEdge(Edge& edge)
{
  SetConstrainedEdge(edge.p, edge.q);
}

// Mark edge as constrained
void Triangle::SetConstrainedEdge(const Point* p, const Point* q)
{
  if ((q == points_[0] && p == points_[1]) || (q == points_[1] && p == points_[0])) {
    constrained_edge_[2] = true;
  } else if ((q == points_[0] && p == points_[2]) || (q == points_[2] && p == points_[0])) {
    constrained_edge_[1] = true;
  } else if ((q == points_[1] && p == points_[2]) || (q == points_[2] && p == points_[1])) {
    constrained_edge_[0] = true;
  }
}

// The point clockwise to given point
const Point* Triangle::PointCW(const Point* point)
{
  if (point == points_[0]) {
    return points_[2];
  } else if (point == points_[1]) {
    return points_[0];
  } else if (point == points_[2]) {
    return points_[1];
  }
  assert(0);
  return nullptr;
}

// The point counter-clockwise to given point
const Point* Triangle::PointCCW(const Point* point)
{
  if (point == points_[0]) {
    return points_[1];
  } else if (point == points_[1]) {
    return points_[2];
  } else if (point == points_[2]) {
    return points_[0];
  }
  assert(0);
  return nullptr;
}

// The neighbor across to given point
Triangle* Triangle::NeighborAcross(const Point* point)
{
  if (point == points_[0]) {
    return neighbors_[0];
  } else if (point == points_[1]) {
    return neighbors_[1];
  } else {
    assert(point == points_[2]);
    return neighbors_[2];
  }
}

// The neighbor clockwise to given point
Triangle* Triangle::NeighborCW(const Point* point)
{
  if (point == points_[0]) {
    return neighbors_[1];
  } else if (point == points_[1]) {
    return neighbors_[2];
  } else {
    assert(point == points_[2]);
    return neighbors_[0];
  }
}

// The neighbor counter-clockwise to given point
Triangle* Triangle::NeighborCCW(const Point* point)
{
  if (point == points_[0]) {
    return neighbors_[2];
  } else if (point == points_[1]) {
    return neighbors_[0];
  } else {
    assert(point == points_[2]);
    return neighbors_[1];
  }
}

bool Triangle::IsConstrainedEdge(int index)
{
  assert(0 <= index && index < 3);
  return constrained_edge_[index];
}

bool Triangle::IsConstrainedEdge(const Point* p)
{
  if (p == points_[0]) {
    return constrained_edge_[0];
  } else if (p == points_[1]) {
    return constrained_edge_[1];
  } else {
    assert(p == points_[2]);
    return constrained_edge_[2];
  }
}

bool Triangle::IsConstrainedEdgeCCW(const Point* p)
{
  if (p == points_[0]) {
    return constrained_edge_[2];
  } else if (p == points_[1]) {
    return constrained_edge_[0];
  } else {
    assert(p == points_[2]);
    return constrained_edge_[1];
  }
}

bool Triangle::IsConstrainedEdgeCW(const Point* p)
{
  if (p == points_[0]) {
    return constrained_edge_[1];
  } else if (p == points_[1]) {
    return constrained_edge_[2];
  } else {
    assert(p == points_[2]);
    return constrained_edge_[0];
  }
}


void Triangle::SetConstrainedEdge(int index, bool ce)
{
  assert(0 <= index && index < 3);
  constrained_edge_[index] = ce;
}

void Triangle::SetConstrainedEdge(const Point* p, bool ce)
{
  if (p == points_[0]) {
    constrained_edge_[0] = ce;
  } else if (p == points_[1]) {
    constrained_edge_[1] = ce;
  } else {
    assert(p == points_[2]);
    constrained_edge_[2] = ce;
  }
}

void Triangle::SetConstrainedEdgeCCW(const Point* p, bool ce)
{
  if (p == points_[0]) {
    constrained_edge_[2] = ce;
  } else if (p == points_[1]) {
    constrained_edge_[0] = ce;
  } else {
    assert(p == points_[2]);
    constrained_edge_[1] = ce;
  }
}

void Triangle::SetConstrainedEdgeCW(const Point* p, bool ce)
{
  if (p == points_[0]) {
    constrained_edge_[1] = ce;
  } else if (p == points_[1]) {
    constrained_edge_[2] = ce;
  } else {
    assert(p == points_[2]);
    constrained_edge_[0] = ce;
  }
}


bool Triangle::IsDelaunayEdge(int index)
{
  assert(0 <= index && index < 3);
  return delaunay_edge_[index];
}

bool Triangle::IsDelaunayEdge(const Point* p)
{
  if (p == points_[0]) {
    return delaunay_edge_[0];
  } else if (p == points_[1]) {
    return delaunay_edge_[1];
  } else {
    assert(p == points_[2]);
    return delaunay_edge_[2];
  }
}

void Triangle::SetDelaunayEdge(int index, bool e)
{
  assert(0 <= index && index < 3);
  assert(!e || neighbors_[index]);
  delaunay_edge_[index] = e;
  if (neighbors_[index]) { neighbors_[index]->SetDelaunayEdge(points_[(index + 1) % 3], points_[(index + 2) % 3], e); }
}

void Triangle::SetDelaunayEdge(const Point* p, bool e)
{
  if (p == points_[0]) {
    assert(!e || neighbors_[0]);
    delaunay_edge_[0] = e;
    if (neighbors_[0]) { neighbors_[0]->SetDelaunayEdge(points_[1], points_[2], e); }
  } else if (p == points_[1]) {
    assert(!e || neighbors_[1]);
    delaunay_edge_[1] = e;
    if (neighbors_[1]) { neighbors_[1]->SetDelaunayEdge(points_[2], points_[0], e); }
  } else {
    assert(p == points_[2]);
    assert(!e || neighbors_[2]);
    delaunay_edge_[2] = e;
    if (neighbors_[2]) { neighbors_[2]->SetDelaunayEdge(points_[0], points_[1], e); }
  }
}

// Does not set the flag of the neighbor triangle: responsability of the caller
void Triangle::SetDelaunayEdge(const Point* p1, const Point* p2, bool e)
{
  if ((p1 == points_[2] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[2]))
    delaunay_edge_[0] = e;
  else if ((p1 == points_[0] && p2 == points_[2]) || (p1 == points_[2] && p2 == points_[0]))
    delaunay_edge_[1] = e;
  else if ((p1 == points_[0] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[0]))
    delaunay_edge_[2] = e;
  else
    assert(0);
}

bool Triangle::CircumcircleContains(const Point& point) const
{
  assert(GetOrientation() == Orientation::CCW);
  const double dx = points_[0]->x - point.x;
  const double dy = points_[0]->y - point.y;
  const double ex = points_[1]->x - point.x;
  const double ey = points_[1]->y - point.y;
  const double fx = points_[2]->x - point.x;
  const double fy = points_[2]->y - point.y;

  const double ap = dx * dx + dy * dy;
  const double bp = ex * ex + ey * ey;
  const double cp = fx * fx + fy * fy;

  return (dx * (fy * bp - cp * ey) - dy * (fx * bp - cp * ex) + ap * (fx * ey - fy * ex)) < 0;
}

Orientation Triangle::GetOrientation() const
{
  return Orient2d(*points_[0], *points_[1], *points_[2]);
}

} // namespace p2t
