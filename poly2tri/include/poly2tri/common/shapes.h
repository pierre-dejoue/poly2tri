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

#include "dll_symbol.h"
#include "point.h"

#include <cassert>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace p2t {

P2T_DLL_SYMBOL std::ostream& operator<<(std::ostream&, const Point&);

// Represents a simple polygon's edge
struct P2T_DLL_SYMBOL Edge {

  const Point *p;   // The lower point
  const Point *q;   // The upper point

  /// Constructor
  Edge(const Point* p1, const Point* p2);
};

P2T_DLL_SYMBOL std::ostream& operator<<(std::ostream&, const Edge&);

enum Orientation { CW, CCW, COLLINEAR };

P2T_DLL_SYMBOL std::ostream& operator<<(std::ostream&, Orientation);

struct Node;

// Triangle-based data structures are know to have better performance than quad-edge structures
// See: J. Shewchuk, "Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator"
//      "Triangulations in CGAL"
class P2T_DLL_SYMBOL Triangle {
public:

  /// Constructor
  Triangle(const Point* a, const Point* b, const Point* c);

  const Point* GetPoint(int index) const;
  const Point* PointCW(const Point* point);
  const Point* PointCCW(const Point* point);
  const Point* OppositePoint(Triangle& t, const Point* p);

  Triangle* GetNeighbor(int index);
  void MarkNeighbor(int index, Triangle& ot);
  void MarkNeighbor(Triangle& ot, int& i, int& oi);
  void ClearNeighbors();

  int Index(const Point* p);
  int EdgeIndex(const Point* p1, const Point* p2);

  Triangle* Neighbor(int index);
  Triangle* NeighborAcross(const Point* point);
  Triangle* NeighborCW(const Point* point);
  Triangle* NeighborCCW(const Point* point);

  bool IsConstrainedEdge(int index);
  bool IsConstrainedEdge(const Point* p);
  bool IsConstrainedEdgeCCW(const Point* p);
  bool IsConstrainedEdgeCW(const Point* p);
  void SetConstrainedEdge(int index, bool ce);
  int SetConstrainedEdge(const Point* p, bool ce);
  int SetConstrainedEdge(Edge& edge);
  int SetConstrainedEdge(const Point* p, const Point* q);

  bool IsDelaunayEdge(int index);
  bool IsDelaunayEdge(const Point* p);
  void SetDelaunayEdge(int index, bool e);
  void SetDelaunayEdge(const Point* p, bool e);

  Node* GetNode(int index);
  Node* GetNode(const Point* p);
  void SetNode(Node& node);
  void ResetNode(Node& node);

  bool Contains(const Point* p);
  bool Contains(const Edge& e);
  bool Contains(const Point* p, const Point* q);

  void Legalize(const Point* p, const Point* op, bool delaunay_edge);
  void Legalize(int index, const Point* op, bool delaunay_edge);

  void ClearDelaunayEdges();

  inline bool IsInterior() const;
  inline void IsInterior(bool b);

  bool CircumcircleContains(const Point&) const;

  Orientation GetOrientation() const;

private:

  void ClearNeighbor(const Triangle* triangle);

  /// Triangle points
  const Point* points_[3];

  /// Backlink to the sweepline nodes that are associated with this triangle
  /// Invariant: node.triangle != nullptr => For some index i, node.p == node.triangle->points_[i] && node.triangle->nodes_[i] == &node
  Node* nodes_[3];

  /// Neighbor list
  Triangle* neighbors_[3];

  /// Flags to determine if an edge is a Constrained edge
  bool constrained_edge_[3];

  /// Flags to determine if an edge is a Delaunay edge
  bool delaunay_edge_[3];

  /// Has this triangle been marked as an interior triangle, meaning a triangle that belongs to the finalized CDT
  bool interior_;

};

P2T_DLL_SYMBOL std::ostream& operator<<(std::ostream&, const Triangle&);

/**
 * Rotates a triangle pair one vertex CW.
 *<pre>
 *
 *  P +-----+ q               P +-----+ q
 *    | t  /|                   |\  t |
 *    |   / |                   | \   |
 *  t1|  /  |t2               t1|  \  |t2
 *    | /   |     after CW:     |   \ |
 *    |/ ot |                   | ot \|
 * oq +-----+ op             oq +-----+ op
 *
 * </pre>
 * Flag delaunay_pair is set to true if the caller is confident the triangle rotation creates a Delaunay pair.
 * Else, it must be set to false, that is for example the case during a FlipScan during an EdgeEvent.
 */
void RotateTrianglePair(Triangle& t, int p, Triangle& ot, int op, bool delaunay_pair = false);

inline const Point* Triangle::GetPoint(int index) const
{
  assert(index < 3);
  return points_[index];
}

inline Triangle* Triangle::GetNeighbor(int index)
{
  assert(index < 3);
  return neighbors_[index];
}

inline bool Triangle::Contains(const Point* p)
{
  return p == points_[0] || p == points_[1] || p == points_[2];
}

inline bool Triangle::Contains(const Edge& e)
{
  return Contains(e.p) && Contains(e.q);
}

inline bool Triangle::Contains(const Point* p, const Point* q)
{
  return Contains(p) && Contains(q);
}

inline bool Triangle::IsInterior() const
{
  return interior_;
}

inline void Triangle::IsInterior(bool b)
{
  interior_ = b;
}

}
