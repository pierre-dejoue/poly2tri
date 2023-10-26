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
/**
 * Sweep-line, Constrained Delaunay Triangulation (CDT) See: Domiter, V. and
 * Zalik, B.(2008)'Sweep-line algorithm for constrained Delaunay triangulation',
 * International Journal of Geographical Information Science
 *
 * "FlipScan" Constrained Edge Algorithm invented by Thomas Åhlén, thahlen@gmail.com
 */

#pragma once

#include <poly2tri/common/shapes.h>
#include <poly2tri/sweep/cdt.h>
#include <poly2tri/sweep/policy.h>

#include <memory>
#include <string_view>
#include <vector>

namespace p2t {

class AdvancingFront;
using BackFront = AdvancingFront;
class SweepContext;
struct Node;
struct Basin;

class Sweep {
public:

  /**
   * Constructor
   */
  Sweep(SweepContext& tcx, CDT::Info& info);

  /**
   * Destructor - clean up memory
   */
  ~Sweep();

  /**
   * Triangulate
   */
  void Triangulate(Policy policy);

private:

  /**
   * Create the advancing front
   */
  AdvancingFront* CreateAdvancingFront();

  /**
   * Create the back front
   *
   * Used to finalize the convex hull triangulation
   */
  BackFront* CreateBackFront();

  /**
   * Delete the advancing ot the back front
   */
  void DeleteFront();

  /**
   * Start sweeping the Y-sorted point set from bottom to top
   */
  void SweepPoints();

  /**
   * Find closes node to the left of the new point and
   * create a new triangle. If needed new holes and basins
   * will be filled to.
   *
   * @param point
   * @return
   */
  Node& PointEvent(const Point* point);

   /**
     * Edge event
     *
     * @param edge
     * @param node
     */
  void EdgeEvent(const Edge* edge, Node* node);

   /**
     * Edge event
     *
     * @param ep
     * @param eq
     * @param triangle
     * @param point
     */
  void EdgeEvent(const Point* ep, const Point* eq, Triangle* triangle, const Point* point);

  /**
   * Creates a new front triangle and legalize it
   *
   * @param point
   * @param node
   * @return
   */
  Node& NewFrontTriangle(const Point* point, Node& node);

  /**
   * Adds a triangle to the sweepline to fill a hole.
   *
   * @param node - IN/OUT. middle node, that is the bottom of the hole.
   *                       It will be removed from the sweepline and
   *                       the node pointer set to *node->prev
   */
  void Fill(Node** node);

  /**
   * Stack the Lawson's legalization of a triangle
   *
   * @param t - The triangle to legalize
   */
  void LegalizePush(Triangle& triangle);

  /**
   * Perform the Lawson's legalization of the triangles on the stack.
   *
   * Consider each neighbor of the triangle and perform the circumcirle check on the pair of triangles.
   * If the check fails, flip both triangles. The process continues recursively until the CDT criteria
   * is valid in the area around the initial triangle.
   */
  void Legalize();

  /**
   * Fills holes in the Advancing Front
   *
   * @param n
   */
  void FillAdvancingFront(Node& n);

  /**
   * Fills a basin that has formed on the Advancing Front to the right
   * of given node.<br>
   * First we decide a left,bottom and right node that forms the
   * boundaries of the basin. Then we fill the basin.
   *
   * @param node - starting node
   */
  void FillRightBasin(Node& node);

  /**
   * Equivalent of FillRightBasin but on the left side of the starting node
   *
   * @param node - starting node
   */
  void FillLeftBasin(Node& node);

  /**
   * Fill a Basin delimited by two monotone chains of vertices
   *
   * @param basin - The description of the basin
   */
  void FillBasin(Basin& basin);

  bool IsEdgeSideOfTriangle(Triangle& triangle, const Point* ep, const Point* eq) const;

  void FillEdgeEvent(const Edge* edge, Node* node);

  void FillRightAboveEdgeEvent(const Edge* edge, Node* node);

  void FillRightBelowEdgeEvent(const Edge* edge, Node& node);

  void FillRightConcaveEdgeEvent(const Edge* edge, Node& node);

  void FillRightConvexEdgeEvent(const Edge* edge, Node& node);

  void FillLeftAboveEdgeEvent(const Edge* edge, Node* node);

  void FillLeftBelowEdgeEvent(const Edge* edge, Node& node);

  void FillLeftConcaveEdgeEvent(const Edge* edge, Node& node);

  void FillLeftConvexEdgeEvent(const Edge* edge, Node& node);

  void FlipEdgeEvent(const Point* ep, const Point* eq, Triangle* t, const Point* p);

   /**
     * When we need to traverse from one triangle to the next we need
     * the point in current triangle that is the opposite point to the next
     * triangle.
     *
     * @param ep
     * @param eq
     * @param ot
     * @param op
     * @return
     */
  const Point* NextFlipPoint(const Point* ep, const Point* eq, Triangle& ot, const Point* op);

   /**
     * Scan part of the FlipScan algorithm<br>
     * When a triangle pair isn't flippable we will scan for the next
     * point that is inside the flip triangle scan area. When found
     * we generate a new flipEdgeEvent
     *
     * @param ep - last point on the edge we are traversing
     * @param eq - first point on the edge we are traversing
     * @param flipTriangle - the current triangle sharing the point eq with edge
     * @param t
     * @param p
     */
  void FlipScanEdgeEvent(const Point* ep, const Point* eq, Triangle& flip_triangle, Triangle& t, const Point* p);

  void FinalizationConvexHull();

  void FinalizationOuterPolygon();

  BackFront* MeshClearBackFrontTriangles(Triangle* tail_triangle);

  void ConvexHullFillOfFront(AdvancingFront& front);    // AdvancingFront or BackFront

  struct EdgeEventData {
    Edge constrained_edge;
    bool right;

    EdgeEventData(const Edge& e, bool r) : constrained_edge(e), right(r) {}
  };

private:
  struct PendingLegalization {
    PendingLegalization(Triangle* t, unsigned int d) : triangle(t), depth(d) {}

    Triangle* triangle;
    unsigned int depth;
  };

  void HandleError(std::string_view msg);

  Node* NewNode(const Point* p, Triangle* t = nullptr);

  // Sweep context
  SweepContext& tcx_;

  // Advancing front. Also used for the back front in the finalization phase.
  std::unique_ptr<AdvancingFront> front_;

  // Nodes of the advancing front
  std::vector<std::unique_ptr<Node>> nodes_;

  // Discarded nodes can be reused
  Node* discarded_nodes_;

  std::unique_ptr<EdgeEventData> edge_event_;

  std::vector<PendingLegalization> legalize_stack_;

  CDT::Info& info_;

};

}
