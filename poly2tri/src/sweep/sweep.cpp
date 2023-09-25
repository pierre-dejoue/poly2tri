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

#include "sweep.h"

#include "advancing_front.h"
#include "back_front.h"
#include "sweep_context.h"
#include "../common/utils.h"

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace p2t {

Sweep::Sweep(SweepContext& tcx) :
  tcx_(tcx),
  front_(),
  nodes_(),
  basin_(),
  edge_event_()
{
}

// Triangulate simple polygon with holes
void Sweep::Triangulate(Policy policy)
{
  CreateAdvancingFront();
  SweepPoints();

  // Finalize the triangulation
  switch (policy)
  {
    case Policy::ConvexHull:
      FinalizationConvexHull();
      break;

    case Policy::OuterPolygon:
    default:
      FinalizationOuterPolygon();
      break;
  }
}

void Sweep::CreateAdvancingFront()
{
  // Initial triangle
  assert(tcx_.points_.size() > 0);
  const Point* lowest_point = tcx_.GetPoint(0);
  Triangle* triangle = tcx_.AddTriangleToMap(lowest_point, tcx_.head(), tcx_.tail());

  Node* af_head = NewNode(tcx_.head(), triangle);
  Node* af_middle = NewNode(lowest_point, triangle);
  Node* af_tail = NewNode(tcx_.tail());

  af_head->next = af_middle;
  af_middle->next = af_tail;
  af_middle->prev = af_head;
  af_tail->prev = af_middle;

  assert(front_ == nullptr);
  front_ = std::make_unique<AdvancingFront>(*af_head, *af_tail);
}

void Sweep::CreateBackFront()
{
  // The back front is oriented backward (from the tail point to the head point)
  Node* bf_head = NewNode(tcx_.tail());
  Node* bf_tail = NewNode(tcx_.head());

  bf_head->next = bf_tail;
  bf_tail->prev = bf_head;

  assert(back_front_ == nullptr);
  back_front_ = std::make_unique<BackFront>(*bf_head, *bf_tail);
}

void Sweep::SweepPoints()
{
  for (size_t i = 1; i < tcx_.point_count(); i++) {
    Node* node = &PointEvent(tcx_.GetPoint(i));
    for (auto* e : tcx_.GetUpperEdges(i)) {
      EdgeEvent(e, node);
    }
  }
}

namespace {

  void FloodFillOfInteriorTriangles(Triangle& interior_triangle)
  {
    std::vector<Triangle *> triangles;
    triangles.push_back(&interior_triangle);

    while(!triangles.empty()) {
      Triangle *t = triangles.back();
      triangles.pop_back();

      if (t != nullptr && !t->IsInterior()) {
        t->IsInterior(true);
        for (int i = 0; i < 3; i++) {
          if (!t->constrained_edge[i])
            triangles.push_back(t->GetNeighbor(i));
        }
      }
    }
  }

}

void Sweep::FinalizationConvexHull()
{
  // Clean the mesh from the two artificial points (head and tail), and simultaneously build the back front
  MeshClearBackFrontTriangles();
  assert(front_); assert(back_front_);

  // Remove exterior triangles from the mesh
  assert(front_->head()->triangle && front_->head()->triangle->IsInterior() == false);
  assert(front_->tail()->triangle == nullptr || front_->tail()->triangle->IsInterior() == false);
  assert(front_->tail()->prev && front_->tail()->prev->triangle && front_->tail()->prev->triangle->IsInterior() == false);
  front_->head()->triangle = nullptr;
  front_->tail()->triangle = nullptr;
  front_->tail()->prev->triangle = nullptr;
  tcx_.MeshCleanExteriorTriangles();

  // Add the bordering triangles to form the convex hull
  ConvexHullFillOfFront(*front_);
  ConvexHullFillOfFront(*back_front_);
}

void Sweep::FinalizationOuterPolygon()
{
  // Get an internal triangle to start with
  Triangle* t = front_->head()->next->triangle;
  const Point* p = front_->head()->next->point;
  while (t && !t->GetConstrainedEdgeCW(p)) {
    t = t->NeighborCCW(p);
  }

  // Collect interior triangles constrained by edges
  if (t)
    FloodFillOfInteriorTriangles(*t);

  // Remove exterior triangles
  tcx_.MeshCleanExteriorTriangles();
}

void Sweep::MeshClearBackFrontTriangles()
{
  CreateBackFront();
  Node* const h = back_front_->head();
  Node* const t = back_front_->tail();
  Node* node = h;
  Node* middle_node = nullptr;

  // Mark "interior" all triangles
  for (auto& tri : tcx_.map_) { tri->IsInterior(true); }

  // Mark "exterior" the triangles of the back front and simultaneously build the back front linked list.
  TraverseBackTriangles(*front_, [this, &node, &middle_node](Triangle* tri, int pivot, bool mid, bool last) {
    tri->IsInterior(false);
    const Point* p = tri->PointCCW(tri->GetPoint(pivot));
    assert(p != back_front_->head()->point && p != back_front_->tail()->point);   // Do not use variables h, t here to prevent a warning in Release
    if (p != node->point) {
      Node* new_node = NewNode(p);
      node->next = new_node;
      new_node->prev = node;
      node = new_node;
    }
    if (mid) {
      middle_node = node;
    } else {
      node->triangle = tri->GetNeighbor(pivot);
      assert(node->triangle == nullptr || node->triangle->IsInterior());
    }
    if (last && !mid) {
      const Point* q = tri->PointCW(tri->GetPoint(pivot));
      assert(q != back_front_->head()->point && q != back_front_->tail()->point);
      Node* new_node = NewNode(q);
      node->next = new_node;
      new_node->prev = node;
      node = new_node;
    }
  });
  node->next = t;
  t->prev = node;
  if (middle_node != nullptr) {
    middle_node->triangle = middle_node->next->triangle;
  }
  h->triangle = h->next->triangle;
}

void Sweep::ConvexHullFillOfFront(AdvancingFront& front)
{
  // This method will fill the nodes from head()->next->next to tail()->prev->prev
  const auto node_range = GetInnerRange(front);
  Node* const begin_node = node_range.first;
  Node* const end_node = node_range.second;
  Node* const min_prev_node = begin_node->prev;
  Node* node = begin_node;
  while (node != end_node)
  {
    assert(node != nullptr);
    if (Orient2d(*node->prev->point, *node->point, *node->next->point) == CCW) {
      Fill(*node);
      if (node->prev != min_prev_node) { node = node->prev; }  // Stay in range begin_node, end_node
      else { node = node->next; }
    }
    else
    {
      node = node->next;
    }
  }
}

Node& Sweep::PointEvent(const Point* point)
{
  Node* node_ptr = front_->LocateNode(point->x);
  if (!node_ptr || !node_ptr->point || !node_ptr->next || !node_ptr->next->point)
  {
    throw std::runtime_error("PointEvent - null node");
  }

  Node& node = *node_ptr;
  Node& new_node = NewFrontTriangle(point, node);

  // Only need to check +epsilon since point never have smaller
  // x value than node due to how we fetch nodes from the front
  if (point->x <= node.point->x + EPSILON) {
    Fill(node);
  }

  //tcx_.AddNode(new_node);

  FillAdvancingFront(new_node);
  return new_node;
}

void Sweep::EdgeEvent(Edge* edge, Node* node)
{
  edge_event_.constrained_edge = edge;
  edge_event_.right = (edge->p->x > edge->q->x);

  if (IsEdgeSideOfTriangle(*node->triangle, edge->p, edge->q)) {
    return;
  }

  // For now we will do all needed filling
  // TODO: integrate with flip process might give some better performance
  //       but for now this avoid the issue with cases that needs both flips and fills
  FillEdgeEvent(edge, node);
  EdgeEvent(edge->p, edge->q, node->triangle, edge->q);
}

void Sweep::EdgeEvent(const Point* ep, const Point* eq, Triangle* triangle, const Point* point)
{
  if (triangle == nullptr) {
    throw std::runtime_error("EdgeEvent - null triangle");
  }
  if (IsEdgeSideOfTriangle(*triangle, ep, eq)) {
    return;
  }

  const Point* p1 = triangle->PointCCW(point);
  const Orientation o1 = Orient2d(*eq, *p1, *ep);
  if (o1 == COLLINEAR) {
    if (triangle->Contains(eq, p1)) {
      triangle->MarkConstrainedEdge(eq, p1);
      // We are modifying the constraint maybe it would be better to
      // not change the given constraint and just keep a variable for the new constraint
      edge_event_.constrained_edge->q = p1;
      triangle = triangle->NeighborAcross(point);
      EdgeEvent(ep, p1, triangle, p1);
    } else {
      throw std::runtime_error("EdgeEvent - collinear points not supported");
    }
    return;
  }

  const Point* p2 = triangle->PointCW(point);
  const Orientation o2 = Orient2d(*eq, *p2, *ep);
  if (o2 == COLLINEAR) {
    if (triangle->Contains(eq, p2)) {
      triangle->MarkConstrainedEdge(eq, p2);
      // We are modifying the constraint maybe it would be better to
      // not change the given constraint and just keep a variable for the new constraint
      edge_event_.constrained_edge->q = p2;
      triangle = triangle->NeighborAcross(point);
      EdgeEvent(ep, p2, triangle, p2);
    } else {
      throw std::runtime_error("EdgeEvent - collinear points not supported");
    }
    return;
  }

  if (o1 == o2) {
    // Need to decide if we are rotating CW or CCW to get to a triangle
    // that will cross edge
    if (o1 == CW) {
      triangle = triangle->NeighborCCW(point);
    } else {
      triangle = triangle->NeighborCW(point);
    }
    EdgeEvent(ep, eq, triangle, point);
  } else {
    // This triangle crosses constraint so lets flippin start!
    assert(triangle);
    FlipEdgeEvent(ep, eq, triangle, point);
  }
}

bool Sweep::IsEdgeSideOfTriangle(Triangle& triangle, const Point* ep, const Point* eq) const
{
  const int index = triangle.EdgeIndex(ep, eq);

  if (index != -1) {
    triangle.MarkConstrainedEdge(index);
    Triangle* t = triangle.GetNeighbor(index);
    if (t) {
      t->MarkConstrainedEdge(ep, eq);
    }
    return true;
  }
  return false;
}

Node& Sweep::NewFrontTriangle(const Point* point, Node& node)
{
  Triangle* triangle = tcx_.AddTriangleToMap(point, node.point, node.next->point);

  triangle->MarkNeighbor(*node.triangle);

  Node* new_node = NewNode(point);

  new_node->next = node.next;
  new_node->prev = &node;
  node.next->prev = new_node;
  node.next = new_node;

  if (!Legalize(*triangle)) {
    MapTriangleToNodes(*triangle);
  }

  return *new_node;
}

void Sweep::Fill(Node& node)
{
  Triangle* triangle = tcx_.AddTriangleToMap(node.prev->point, node.point, node.next->point);

  // TODO: should copy the constrained_edge value from neighbor triangles
  //       for now constrained_edge values are copied during the legalize
  triangle->MarkNeighbor(*node.prev->triangle);
  triangle->MarkNeighbor(*node.triangle);

  // Update the advancing front
  node.prev->next = node.next;
  node.next->prev = node.prev;

  // If it was legalized the triangle has already been mapped
  if (!Legalize(*triangle)) {
    MapTriangleToNodes(*triangle);
  }
}

void Sweep::FillAdvancingFront(Node& n)
{

  // Fill right holes
  Node* node = n.next;

  while (node && node->next) {
    // if HoleAngle exceeds 90 degrees then break.
    if (LargeHole_DontFill(node)) break;
    Fill(*node);
    node = node->next;
  }

  // Fill left holes
  node = n.prev;

  while (node && node->prev) {
    // if HoleAngle exceeds 90 degrees then break.
    if (LargeHole_DontFill(node)) break;
    Fill(*node);
    node = node->prev;
  }

  // Fill right basins
  if (n.next && n.next->next) {
    const double angle = BasinAngle(n);
    if (angle < PI_3div4) {
      FillBasin(n);
    }
  }
}

// True if HoleAngle exceeds 90 degrees.
// LargeHole_DontFill checks if the advancing front has a large hole.
// A "Large hole" is a triangle formed by a sequence of points in the advancing
// front where three neighbor points form a triangle.
// And angle between left-top, bottom, and right-top points is more than 90 degrees.
// The first part of the algorithm reviews only three neighbor points, e.g. named A, B, C.
// Additional part of this logic reviews a sequence of 5 points -
// additionally reviews one point before and one after the sequence of three (A, B, C),
// e.g. named X and Y.
// In this case, angles are XBC and ABY and this if angles are negative or more
// than 90 degrees LargeHole_DontFill returns true.
// But there is a configuration when ABC has a negative angle but XBC or ABY is less
// than 90 degrees and positive.
// Then function LargeHole_DontFill return false and initiates filling.
// This filling creates a triangle ABC and adds it to the advancing front.
// But in the case when angle ABC is negative this triangle goes inside the advancing front
// and can intersect previously created triangles.
// This triangle leads to making wrong advancing front and problems in triangulation in the future.
// Looks like such a triangle should not be created.
// The simplest way to check and fix it is to check an angle ABC.
// If it is negative LargeHole_DontFill should return true and
// not initiate creating the ABC triangle in the advancing front.
// X______A         Y
//        \        /
//         \      /
//          \ B  /
//           |  /
//           | /
//           |/
//           C
bool Sweep::LargeHole_DontFill(const Node* node)
{
  const Node* nextNode = node->next;
  const Node* prevNode = node->prev;
  if (!AngleExceeds90Degrees(node->point, nextNode->point, prevNode->point))
          return false;

  if (AngleIsNegative(node->point, nextNode->point, prevNode->point))
          return true;

  // Check additional points on front.
  const Node* next2Node = nextNode->next;
  // "..Plus.." because only want angles on same side as point being added.
  if ((next2Node != nullptr) && !AngleExceedsPlus90DegreesOrIsNegative(node->point, next2Node->point, prevNode->point))
          return false;

  const Node* prev2Node = prevNode->prev;
  // "..Plus.." because only want angles on same side as point being added.
  if ((prev2Node != nullptr) && !AngleExceedsPlus90DegreesOrIsNegative(node->point, nextNode->point, prev2Node->point))
          return false;

  return true;
}

bool Sweep::AngleIsNegative(const Point* origin, const Point* pa, const Point* pb)
{
    const double angle = Angle(origin, pa, pb);
    return angle < 0;
}

bool Sweep::AngleExceeds90Degrees(const Point* origin, const Point* pa, const Point* pb)
{
  const double angle = Angle(origin, pa, pb);
  return ((angle > PI_div2) || (angle < -PI_div2));
}

bool Sweep::AngleExceedsPlus90DegreesOrIsNegative(const Point* origin, const Point* pa, const Point* pb)
{
  const double angle = Angle(origin, pa, pb);
  return (angle > PI_div2) || (angle < 0);
}

double Sweep::Angle(const Point* origin, const Point* pa, const Point* pb)
{
  /* Complex plane
   * ab = cosA +i*sinA
   * ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
   * atan2(y,x) computes the principal value of the argument function
   * applied to the complex number x+iy
   * Where x = ax*bx + ay*by
   *       y = ax*by - ay*bx
   */
  const double px = origin->x;
  const double py = origin->y;
  const double ax = pa->x - px;
  const double ay = pa->y - py;
  const double bx = pb->x - px;
  const double by = pb->y - py;
  const double x = ax * by - ay * bx;
  const double y = ax * bx + ay * by;
  return atan2(x, y);
}

double Sweep::BasinAngle(const Node& node)
{
  const double ax = node.point->x - node.next->next->point->x;
  const double ay = node.point->y - node.next->next->point->y;
  return atan2(ay, ax);
}

double Sweep::HoleAngle(const Node& node)
{
  /* Complex plane
   * ab = cosA +i*sinA
   * ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
   * atan2(y,x) computes the principal value of the argument function
   * applied to the complex number x+iy
   * Where x = ax*bx + ay*by
   *       y = ax*by - ay*bx
   */
  const double ax = node.next->point->x - node.point->x;
  const double ay = node.next->point->y - node.point->y;
  const double bx = node.prev->point->x - node.point->x;
  const double by = node.prev->point->y - node.point->y;
  return atan2(ax * by - ay * bx, ax * bx + ay * by);
}

bool Sweep::Legalize(Triangle& t)
{
  // To legalize a triangle we start by finding if any of the three edges
  // violate the Delaunay condition
  for (int i = 0; i < 3; i++) {
    if (t.delaunay_edge[i])
      continue;

    Triangle* ot = t.GetNeighbor(i);

    if (ot) {
      const Point* p = t.GetPoint(i);
      const Point* op = ot->OppositePoint(t, p);
      int oi = ot->Index(op);

      // If this is a Constrained Edge or a Delaunay Edge(only during recursive legalization)
      // then we should not try to legalize
      if (ot->constrained_edge[oi] || ot->delaunay_edge[oi]) {
        t.constrained_edge[i] = ot->constrained_edge[oi];
        continue;
      }

      bool inside = Incircle(*p, *t.PointCCW(p), *t.PointCW(p), *op);

      if (inside) {
        // Lets mark this shared edge as Delaunay
        t.delaunay_edge[i] = true;
        ot->delaunay_edge[oi] = true;

        // Lets rotate shared edge one vertex CW to legalize it
        RotateTrianglePair(t, p, *ot, op);

        // We now got one valid Delaunay Edge shared by two triangles
        // This gives us 4 new edges to check for Delaunay

        // Make sure that triangle to node mapping is done only one time for a specific triangle
        bool not_legalized = !Legalize(t);
        if (not_legalized) {
          MapTriangleToNodes(t);
        }

        not_legalized = !Legalize(*ot);
        if (not_legalized)
          MapTriangleToNodes(*ot);

        // Reset the Delaunay edges, since they only are valid Delaunay edges
        // until we add a new triangle or point.
        // XXX: need to think about this. Can these edges be tried after we
        //      return to previous recursive level?
        t.delaunay_edge[i] = false;
        ot->delaunay_edge[oi] = false;

        // If triangle have been legalized no need to check the other edges since
        // the recursive legalization will handles those so we can end here.
        return true;
      }
    }
  }
  return false;
}

bool Sweep::Incircle(const Point& pa, const Point& pb, const Point& pc, const Point& pd)
{
  const double adx = pa.x - pd.x;
  const double ady = pa.y - pd.y;
  const double bdx = pb.x - pd.x;
  const double bdy = pb.y - pd.y;

  const double adxbdy = adx * bdy;
  const double bdxady = bdx * ady;
  const double oabd = adxbdy - bdxady;

  if (oabd <= 0)
    return false;

  const double cdx = pc.x - pd.x;
  const double cdy = pc.y - pd.y;

  const double cdxady = cdx * ady;
  const double adxcdy = adx * cdy;
  const double ocad = cdxady - adxcdy;

  if (ocad <= 0)
    return false;

  const double bdxcdy = bdx * cdy;
  const double cdxbdy = cdx * bdy;
  const double obcd = bdxcdy - cdxbdy;

  const double alift = adx * adx + ady * ady;
  const double blift = bdx * bdx + bdy * bdy;
  const double clift = cdx * cdx + cdy * cdy;

  const double det = alift * obcd + blift * ocad + clift * oabd;

  return det > 0;
}

void Sweep::RotateTrianglePair(Triangle& t, const Point* p, Triangle& ot, const Point* op)
{
  Triangle* n1, *n2, *n3, *n4;
  n1 = t.NeighborCCW(p);
  n2 = t.NeighborCW(p);
  n3 = ot.NeighborCCW(op);
  n4 = ot.NeighborCW(op);

  bool ce1, ce2, ce3, ce4;
  ce1 = t.GetConstrainedEdgeCCW(p);
  ce2 = t.GetConstrainedEdgeCW(p);
  ce3 = ot.GetConstrainedEdgeCCW(op);
  ce4 = ot.GetConstrainedEdgeCW(op);

  bool de1, de2, de3, de4;
  de1 = t.GetDelaunayEdgeCCW(p);
  de2 = t.GetDelaunayEdgeCW(p);
  de3 = ot.GetDelaunayEdgeCCW(op);
  de4 = ot.GetDelaunayEdgeCW(op);

  t.Legalize(p, op);
  ot.Legalize(op, p);

  // Remap delaunay_edge
  ot.SetDelaunayEdgeCCW(p, de1);
  t.SetDelaunayEdgeCW(p, de2);
  t.SetDelaunayEdgeCCW(op, de3);
  ot.SetDelaunayEdgeCW(op, de4);

  // Remap constrained_edge
  ot.SetConstrainedEdgeCCW(p, ce1);
  t.SetConstrainedEdgeCW(p, ce2);
  t.SetConstrainedEdgeCCW(op, ce3);
  ot.SetConstrainedEdgeCW(op, ce4);

  // Remap neighbors
  // XXX: might optimize the markNeighbor by keeping track of
  //      what side should be assigned to what neighbor after the
  //      rotation. Now mark neighbor does lots of testing to find
  //      the right side.
  t.ClearNeighbors();
  ot.ClearNeighbors();
  if (n1) ot.MarkNeighbor(*n1);
  if (n2) t.MarkNeighbor(*n2);
  if (n3) t.MarkNeighbor(*n3);
  if (n4) ot.MarkNeighbor(*n4);
  t.MarkNeighbor(ot);
}

void Sweep::FillBasin(Node& node)
{
  if (Orient2d(*node.point, *node.next->point, *node.next->next->point) == CCW) {
    basin_.left_node = node.next->next;
  } else {
    basin_.left_node = node.next;
  }

  // Find the bottom and right node
  basin_.bottom_node = basin_.left_node;
  while (basin_.bottom_node->next
         && basin_.bottom_node->point->y >= basin_.bottom_node->next->point->y) {
    basin_.bottom_node = basin_.bottom_node->next;
  }
  if (basin_.bottom_node == basin_.left_node) {
    // No valid basin
    return;
  }

  basin_.right_node = basin_.bottom_node;
  while (basin_.right_node->next
         && basin_.right_node->point->y < basin_.right_node->next->point->y) {
    basin_.right_node = basin_.right_node->next;
  }
  if (basin_.right_node == basin_.bottom_node) {
    // No valid basins
    return;
  }

  basin_.width = basin_.right_node->point->x - basin_.left_node->point->x;
  basin_.left_highest = basin_.left_node->point->y > basin_.right_node->point->y;

  FillBasinReq(basin_.bottom_node);
}

void Sweep::FillBasinReq(Node* node)
{
  // if shallow stop filling
  if (IsShallow(*node)) {
    return;
  }

  Fill(*node);

  if (node->prev == basin_.left_node && node->next == basin_.right_node) {
    return;
  } else if (node->prev == basin_.left_node) {
    Orientation o = Orient2d(*node->point, *node->next->point, *node->next->next->point);
    if (o == CW) {
      return;
    }
    node = node->next;
  } else if (node->next == basin_.right_node) {
    Orientation o = Orient2d(*node->point, *node->prev->point, *node->prev->prev->point);
    if (o == CCW) {
      return;
    }
    node = node->prev;
  } else {
    // Continue with the neighbor node with lowest Y value
    if (node->prev->point->y < node->next->point->y) {
      node = node->prev;
    } else {
      node = node->next;
    }
  }

  FillBasinReq(node);
}

bool Sweep::IsShallow(Node& node) const
{
  double height;

  if (basin_.left_highest) {
    height = basin_.left_node->point->y - node.point->y;
  } else {
    height = basin_.right_node->point->y - node.point->y;
  }

  // if shallow stop filling
  return basin_.width > height;
}

void Sweep::FillEdgeEvent(Edge* edge, Node* node)
{
  if (edge_event_.right) {
    FillRightAboveEdgeEvent(edge, node);
  } else {
    FillLeftAboveEdgeEvent(edge, node);
  }
}

void Sweep::FillRightAboveEdgeEvent(Edge* edge, Node* node)
{
  while (node->next->point->x < edge->p->x) {
    // Check if next node is below the edge
    if (Orient2d(*edge->q, *node->next->point, *edge->p) == CCW) {
      FillRightBelowEdgeEvent(edge, *node);
    } else {
      node = node->next;
    }
  }
}

void Sweep::FillRightBelowEdgeEvent(Edge* edge, Node& node)
{
  if (node.point->x < edge->p->x) {
    if (Orient2d(*node.point, *node.next->point, *node.next->next->point) == CCW) {
      // Concave
      FillRightConcaveEdgeEvent(edge, node);
    } else {
      // Convex
      FillRightConvexEdgeEvent(edge, node);
      // Retry this one
      FillRightBelowEdgeEvent(edge, node);
    }
  }
}

void Sweep::FillRightConcaveEdgeEvent(Edge* edge, Node& node)
{
  Fill(*node.next);
  if (node.next->point != edge->p) {
    // Next above or below edge?
    if (Orient2d(*edge->q, *node.next->point, *edge->p) == CCW) {
      // Below
      if (Orient2d(*node.point, *node.next->point, *node.next->next->point) == CCW) {
        // Next is concave
        FillRightConcaveEdgeEvent(edge, node);
      } else {
        // Next is convex
      }
    }
  }
}

void Sweep::FillRightConvexEdgeEvent(Edge* edge, Node& node)
{
  // Next concave or convex?
  if (Orient2d(*node.next->point, *node.next->next->point, *node.next->next->next->point) == CCW) {
    // Concave
    FillRightConcaveEdgeEvent(edge, *node.next);
  } else {
    // Convex
    // Next above or below edge?
    if (Orient2d(*edge->q, *node.next->next->point, *edge->p) == CCW) {
      // Below
      FillRightConvexEdgeEvent(edge, *node.next);
    } else {
      // Above
    }
  }
}

void Sweep::FillLeftAboveEdgeEvent(Edge* edge, Node* node)
{
  while (node->prev->point->x > edge->p->x) {
    // Check if next node is below the edge
    if (Orient2d(*edge->q, *node->prev->point, *edge->p) == CW) {
      FillLeftBelowEdgeEvent(edge, *node);
    } else {
      node = node->prev;
    }
  }
}

void Sweep::FillLeftBelowEdgeEvent(Edge* edge, Node& node)
{
  if (node.point->x > edge->p->x) {
    if (Orient2d(*node.point, *node.prev->point, *node.prev->prev->point) == CW) {
      // Concave
      FillLeftConcaveEdgeEvent(edge, node);
    } else {
      // Convex
      FillLeftConvexEdgeEvent(edge, node);
      // Retry this one
      FillLeftBelowEdgeEvent(edge, node);
    }
  }
}

void Sweep::FillLeftConvexEdgeEvent(Edge* edge, Node& node)
{
  // Next concave or convex?
  if (Orient2d(*node.prev->point, *node.prev->prev->point, *node.prev->prev->prev->point) == CW) {
    // Concave
    FillLeftConcaveEdgeEvent(edge, *node.prev);
  } else {
    // Convex
    // Next above or below edge?
    if (Orient2d(*edge->q, *node.prev->prev->point, *edge->p) == CW) {
      // Below
      FillLeftConvexEdgeEvent(edge, *node.prev);
    } else {
      // Above
    }
  }
}

void Sweep::FillLeftConcaveEdgeEvent(Edge* edge, Node& node)
{
  Fill(*node.prev);
  if (node.prev->point != edge->p) {
    // Next above or below edge?
    if (Orient2d(*edge->q, *node.prev->point, *edge->p) == CW) {
      // Below
      if (Orient2d(*node.point, *node.prev->point, *node.prev->prev->point) == CW) {
        // Next is concave
        FillLeftConcaveEdgeEvent(edge, node);
      } else {
        // Next is convex
      }
    }
  }
}

void Sweep::FlipEdgeEvent(const Point* ep, const Point* eq, Triangle* t, const Point* p)
{
  assert(t);
  Triangle* ot_ptr = t->NeighborAcross(p);
  if (ot_ptr == nullptr)
  {
    throw std::runtime_error("FlipEdgeEvent - null neighbor across");
  }
  Triangle& ot = *ot_ptr;
  const Point* op = ot.OppositePoint(*t, p);

  if (InScanArea(*p, *t->PointCCW(p), *t->PointCW(p), *op)) {
    // Lets rotate shared edge one vertex CW
    RotateTrianglePair(*t, p, ot, op);
    MapTriangleToNodes(*t);
    MapTriangleToNodes(ot);

    if (p == eq && op == ep) {
      if (eq == edge_event_.constrained_edge->q && ep == edge_event_.constrained_edge->p) {
        t->MarkConstrainedEdge(ep, eq);
        ot.MarkConstrainedEdge(ep, eq);
        Legalize(*t);
        Legalize(ot);
      } else {
        // XXX: I think one of the triangles should be legalized here?
      }
    } else {
      const Orientation o = Orient2d(*eq, *op, *ep);
      t = &NextFlipTriangle(o, *t, ot, p, op);
      FlipEdgeEvent(ep, eq, t, p);
    }
  } else {
    const Point* new_p = NextFlipPoint(ep, eq, ot, op);
    FlipScanEdgeEvent(ep, eq, *t, ot, new_p);
    EdgeEvent(ep, eq, t, p);
  }
}

Triangle& Sweep::NextFlipTriangle(Orientation o, Triangle& t, Triangle& ot, const Point* p, const Point* op)
{
  if (o == CCW) {
    // ot is not crossing edge after flip
    int edge_index = ot.EdgeIndex(p, op);
    ot.delaunay_edge[edge_index] = true;
    Legalize(ot);
    ot.ClearDelaunayEdges();
    return t;
  }

  // t is not crossing edge after flip
  int edge_index = t.EdgeIndex(p, op);

  t.delaunay_edge[edge_index] = true;
  Legalize(t);
  t.ClearDelaunayEdges();
  return ot;
}

const Point* Sweep::NextFlipPoint(const Point* ep, const Point* eq, Triangle& ot, const Point* op)
{
  const Orientation o2d = Orient2d(*eq, *op, *ep);
  if (o2d == CW) {
    // Right
    return ot.PointCCW(op);
  } else if (o2d == CCW) {
    // Left
    return ot.PointCW(op);
  } else {
    assert(o2d == COLLINEAR);
    throw std::runtime_error("[Unsupported] Opposing point on constrained edge");
  }
}

void Sweep::FlipScanEdgeEvent(const Point* ep, const Point* eq, Triangle& flip_triangle,
                              Triangle& t, const Point* p)
{
  Triangle* ot_ptr = t.NeighborAcross(p);
  if (ot_ptr == nullptr) {
    throw std::runtime_error("FlipScanEdgeEvent - null neighbor across");
  }

  const Point* op = ot_ptr->OppositePoint(t, p);
  if (op == nullptr) {
    throw std::runtime_error("FlipScanEdgeEvent - null opposing point");
  }

  const Point* p1 = flip_triangle.PointCCW(eq);
  const Point* p2 = flip_triangle.PointCW(eq);
  if (p1 == nullptr || p2 == nullptr) {
    throw std::runtime_error("FlipScanEdgeEvent - null on either of points");
  }

  Triangle& ot = *ot_ptr;

  if (InScanArea(*eq, *p1, *p2, *op)) {
    // flip with new edge op->eq
    FlipEdgeEvent(eq, op, &ot, op);
    // TODO: Actually I just figured out that it should be possible to
    //       improve this by getting the next ot and op before the the above
    //       flip and continue the flipScanEdgeEvent here
    // set new ot and op here and loop back to inScanArea test
    // also need to set a new flip_triangle first
    // Turns out at first glance that this is somewhat complicated
    // so it will have to wait.
  } else {
    const Point* new_p = NextFlipPoint(ep, eq, ot, op);
    FlipScanEdgeEvent(ep, eq, flip_triangle, ot, new_p);
  }
}

void Sweep::MapTriangleToNodes(Triangle& t)
{
  for (int i = 0; i < 3; i++) {
    if (!t.GetNeighbor(i)) {
      Node* n = front_->LocatePoint(t.PointCW(t.GetPoint(i)));
      if (n)
        n->triangle = &t;
    }
  }
}

Node* Sweep::NewNode(const Point* p, Triangle* t)
{
  nodes_.emplace_back(std::make_unique<Node>(p, t));
  return nodes_.back().get();
}

Sweep::~Sweep() = default;

} // namespace p2t
