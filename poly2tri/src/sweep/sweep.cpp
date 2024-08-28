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

#include "sweep.h"

#include "advancing_front.h"
#include "back_triangles.h"
#include "sweep_context.h"
#include "trace.h"
#include "../common/node.h"
#include "../common/predicates.h"
#include "../common/utils.h"

#include <poly2tri/common/shapes_io.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <stdexcept>

namespace p2t {

class NodeStorage {
public:
  NodeStorage()
    : capacity{0}
    , next_batch_sz{FIRST_BATCH_SZ}
    , next_batch_idx{0}
    , available_nodes{nullptr}
    , node_buffers{}
  {
    // Memory is allocated in successive batches to ensure the next/prev pointers of the Nodes remain valid.
    // The size of the batches follow a geometric progression starting with a first allocation of FIRST_BATCH_SZ (256) Nodes.
    //
    //      Batch_sz   Total_capacity      Mem_footprint
    //      256        256                 10 kB
    //      256        512                 20 kB
    //      512        1024                40 kB
    //      1024       2048                80 kB
    //      ...        ...                 ...
    //
    PrepareNextBatchOfNodes(FIRST_BATCH_SZ);
  }

  Node* NewNode()
  {
    if (!available_nodes) {
      PrepareNextBatchOfNodes(next_batch_sz);
      next_batch_sz = 2 * next_batch_sz;
    }
    assert(available_nodes);
    Node* new_node = available_nodes;
    available_nodes = new_node->next;
    new_node->next = nullptr;
    new_node->prev = nullptr;
    return new_node;
  }

  void DiscardNode(Node* node)
  {
    assert(node);
    node->next = available_nodes;
    available_nodes = node;
  }

  // Storage total capacity in number of Nodes
  std::size_t Capacity() const
  {
    return capacity;
  }

  std::size_t MemoryFootprint() const
  {
    return capacity * sizeof(Node);
  }

private:
  static constexpr std::size_t FIRST_BATCH_SZ = 256u;

  void PrepareNextBatchOfNodes(std::size_t batch_sz)
  {
    assert(available_nodes == nullptr);
    if (next_batch_idx >= node_buffers.size())
    {
       throw std::runtime_error("NodeStorage: Out of memory capacity");
    }
    assert(!node_buffers[next_batch_idx]);
    node_buffers[next_batch_idx] = std::unique_ptr<Node[]>(new Node[batch_sz]);
    Node* new_node = node_buffers[next_batch_idx++].get();
    std::size_t count = batch_sz;
    Node** node_pptr = &available_nodes;
    while (count--) {
      *node_pptr = new_node;
      node_pptr = &(new_node->next);
      new_node++;
    }
    capacity += batch_sz;
    assert(available_nodes != nullptr);
  }

private:
  std::size_t capacity;
  std::size_t next_batch_sz;
  std::size_t next_batch_idx;

  // Linked list of available nodes. Discarded nodes are put back in the list and can be reused.
  Node* available_nodes;

  std::array<std::unique_ptr<Node[]>, 32> node_buffers;
};


Sweep::Sweep(SweepContext& tcx, CDT::Info& info) :
  tcx_(tcx),
  front_(),
  node_storage_(),
  edge_event_(),
  legalize_stack_(),
  info_(info)
{
  node_storage_ = std::make_unique<NodeStorage>();
}

Sweep::~Sweep() {
  if (node_storage_)
    info_.nodes_memory_footprint_in_bytes = node_storage_->MemoryFootprint();
}

// Triangulate simple polygon with holes
void Sweep::Triangulate(Policy policy)
{
  CreateAdvancingFront();
  SweepPoints();

  info_.nb_triangles_pre_finalization = static_cast<unsigned int>(tcx_.TriangleStorageNbOfCreatedTriangles());

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

  info_.nb_output_triangles = static_cast<unsigned int>(tcx_.GetTriangles().size());
}

AdvancingFront* Sweep::CreateAdvancingFront()
{
  // Initial triangle
  assert(tcx_.points_.size() > 0);
  const Point* lowest_point = tcx_.GetPoint(0);
  Triangle* triangle = tcx_.AddTriangle(lowest_point, tcx_.head(), tcx_.tail());

  Node* af_head = NewNode(tcx_.head(), triangle);
  Node* af_middle = NewNode(lowest_point, triangle);
  Node* af_tail = NewNode(tcx_.tail());

  assert(af_head->prev == nullptr);
  af_head->next = af_middle;
  af_middle->prev = af_head;
  af_middle->next = af_tail;
  af_tail->prev = af_middle;
  assert(af_tail->next == nullptr);

  assert(!front_);
  front_ = std::make_optional<AdvancingFront>(*af_head, *af_tail);

  TRACE_OUT << "CreateAdvancingFront - advancing_front=" << front_.value() << std::endl;

  return &front_.value();
}

BackFront* Sweep::CreateBackFront()
{
  TRACE_OUT << "CreateBackFront" << std::endl;

  // The back front is oriented backward (from the tail point to the head point)
  Node* bf_head = NewNode(tcx_.tail());
  Node* bf_tail = NewNode(tcx_.head());

  assert(bf_head->prev == nullptr);
  bf_head->next = bf_tail;
  bf_tail->prev = bf_head;
  assert(bf_tail->next == nullptr);

  assert(!front_);
  front_ = std::make_optional<BackFront>(*bf_head, *bf_tail);

  return &front_.value();
}

void Sweep::DeleteFront()
{
  assert(front_);
  Node* node = front_->head();
  while (node != nullptr) {
    node->ResetTriangle();
    node->prev = nullptr;
    Node* delete_me = node;
    node = node->next;
    node_storage_->DiscardNode(delete_me);
  }
  front_.reset();
}

void Sweep::SweepPoints()
{
  const std::size_t point_count = tcx_.points_.size();
  // The point with index 0 is already in the advancing front
  for (size_t i = 1; i < point_count; i++) {
    Node* node = &PointEvent(tcx_.GetPoint(i));
    Legalize();
    for (const auto& e : tcx_.GetUpperEdges(i)) {
      EdgeEvent(&e, node);
      Legalize();
    }
  }
}

namespace {

  void FloodFillOfInteriorTriangles(Triangle& interior_triangle)
  {
    std::vector<Triangle*> triangle_stack;
    triangle_stack.emplace_back(&interior_triangle);

    while(!triangle_stack.empty()) {
      Triangle* t = triangle_stack.back();
      triangle_stack.pop_back();
      assert(t);
      if (t->GetState() == Triangle::State::Interior) {
        continue;
      }
      t->SetInterior();
      for (int i = 0; i < 3; i++) {
        if (!t->IsConstrainedEdge(i)) {
          auto* nt = t->GetNeighbor(i);
          if (nt && nt->GetState() == Triangle::State::Normal) {
            triangle_stack.emplace_back(nt);
          }
        }
      }
    }
  }

  bool SetConstrainedEdgeIfSideOfTriangle(Triangle& triangle, const Point* ep, const Point* eq)
  {
    assert(ep); assert(eq);
    const int index = triangle.EdgeIndex(ep, eq);
    const bool is_side_of_triangle = index != -1;
    if (is_side_of_triangle) {
      triangle.SetConstrainedEdge(index, true);
      Triangle* t = triangle.GetNeighbor(index);
      if (t) { t->SetConstrainedEdge(ep, eq); }
    }
    return is_side_of_triangle;
  }

} // namespace

void Sweep::FinalizationConvexHull()
{
  // Add the bordering triangles to form the convex hull of the advancing front
  assert(front_);
  TRACE_OUT << "FinalizationConvexHull -\n"
            << "  advancing_front=" << *front_ << std::endl;
  ConvexHullFillOfFront(*front_);

  // Tail triangle
  assert(front_->tail()->prev);
  Triangle* tail_triangle = front_->tail()->prev->triangle;

  // Delete the advancing front and release the nodes
  DeleteFront();

  // Clean the mesh from the two artificial points (head and tail), and simultaneously build the back front
  auto* back_front = MeshClearBackFrontTriangles(tail_triangle);
  assert(back_front);

  // Add the bordering triangles to form the convex hull of the back front
  TRACE_OUT << "FinalizationConvexHull -\n"
            << "  back_front=" << *back_front << std::endl;
  ConvexHullFillOfFront(*back_front);
  DeleteFront();

  // The final triangulation is made of all the non-discarded Triangles
  tcx_.PopulateTriangleMap(Triangle::State::Normal);
}

void Sweep::FinalizationOuterPolygon()
{
  // Get an internal triangle to start with
  Triangle* t = front_->head()->next->triangle;
  const Point* p = front_->head()->next->point;
  while (t && !t->IsConstrainedEdgeCW(p)) {
    t = t->NeighborCCW(p);
  }

  // Collect interior triangles constrained by edges
  if (t) { FloodFillOfInteriorTriangles(*t); }

  // Delete the advancing front and release the nodes
  DeleteFront();

  // Filter out the external Triangle to get the triangulation
  tcx_.PopulateTriangleMap(Triangle::State::Interior);
}

BackFront* Sweep::MeshClearBackFrontTriangles(Triangle* tail_triangle)
{
  assert(tail_triangle);
  auto* back_front = CreateBackFront();
  Node* const hd = back_front->head();
  Node* const tl = back_front->tail();
  Node* node = hd;

  Node* triangles_to_discard = nullptr;

  // Mark "to_discard" the triangles attached with the head and tail artificial points and simultaneously build the back front linked list.
  TraverseBackTriangles(tcx_.head(), tcx_.tail(), tail_triangle, [this, &node, &triangles_to_discard](Triangle* tri, int pivot, bool mid, bool last) {
    // Mark current triangle "to be discarded"
    Node* discard_node = node_storage_->NewNode();
    discard_node->triangle = tri;
    discard_node->next = triangles_to_discard;
    triangles_to_discard = discard_node;
    // Build the back front
    if (mid && !last)
      return;
    const Point* p = tri->PointCCW(tri->GetPoint(pivot));
    assert(p != front_->head()->point && p != front_->tail()->point);   // Do not use variables hd and tl here to prevent a warning in Release
    if (p != node->point) {
      Node* new_node = NewNode(p);
      node->next = new_node;
      new_node->prev = node;
      node = new_node;
    }
    if (!mid) {
      node->SetTriangle(tri->GetNeighbor(pivot));
    }
    if (last && !mid) {
      const Point* q = tri->PointCW(tri->GetPoint(pivot));
      assert(q != front_->head()->point && q != front_->tail()->point);
      Node* new_node = NewNode(q);
      node->next = new_node;
      new_node->prev = node;
      node = new_node;
    }
  });
  node->next = tl;
  tl->prev = node;

  // Discard the back triangles
  while (triangles_to_discard) {
    assert(triangles_to_discard->triangle);
    tcx_.DiscardTriangle(*triangles_to_discard->triangle);
    triangles_to_discard->triangle = nullptr;
    Node* delete_me = triangles_to_discard;
    triangles_to_discard = triangles_to_discard->next;
    node_storage_->DiscardNode(delete_me);
  }

  // Ensure that all the node triangles in the back front are normal ones
  node = hd;
  while (node != nullptr && node != tl)
  {
    if (node->triangle && node->triangle->GetState() == Triangle::State::Discarded) {
      node->ResetTriangle();
    }
    node = node->next;
  }

  return back_front;
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
      Fill(&node);   // node is set to node->prev
      // Move backward in the list but stay in range:  begin_node <= node < end_node
      if (node == min_prev_node) { node = node->next; }
    }
    else
    {
      node = node->next;
    }
  }
  Legalize();
}

Node& Sweep::PointEvent(const Point* point)
{
  TRACE_OUT << "PointEvent - point=" << *point << std::endl;
  Node* node_ptr = front_->LocateNode(point->x);
  if (!node_ptr || !node_ptr->point || !node_ptr->next || !node_ptr->next->point)
  {
    HandleError("PointEvent - null node");
  }
  assert(node_ptr);
  assert(node_ptr->point);

  Node& node = *node_ptr;
  Node& new_node = NewFrontTriangle(point, node);

  // In case the projection coincides with the node's point, fill the previous node
  // Note: Legitimate use of the equal operator for comparing the input points' coordinates
  if (point->x == node.point->x) {
    TRACE_OUT << "PointEvent - The projection coincides with front node=" << node << std::endl;
    Fill(&node_ptr);
  }

  FillAdvancingFront(new_node);
  return new_node;
}

void Sweep::EdgeEvent(const Edge* edge, Node* node)
{
  assert(edge); assert(node);
  TRACE_OUT << "EdgeEvent - "
            << "edge=" << *edge << "; "
            << "node=" << *node
            << std::endl;

  // TODO: Reset the edge event when done
  edge_event_ = std::make_optional<EdgeEventData>(*edge, edge->p->x > edge->q->x);

  if (SetConstrainedEdgeIfSideOfTriangle(*node->triangle, edge->p, edge->q)) {
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
    HandleError("EdgeEvent - null triangle");
  }
  assert(triangle);
  assert(edge_event_);
  TRACE_OUT << "EdgeEvent - "
            << "edge={ ep=" << *ep << ", eq=" << *eq << " }; "
            << "triangle=" << *triangle << "; "
            << "point=" << *point
            << std::endl;
  if (SetConstrainedEdgeIfSideOfTriangle(*triangle, ep, eq)) {
    return;
  }

  const Point* p1 = triangle->PointCCW(point);
  const Orientation o1 = Orient2d(*eq, *p1, *ep);
  if (o1 == COLLINEAR) {
    TRACE_OUT << "EdgeEvent - o1=" << o1 << std::endl;
    if (SetConstrainedEdgeIfSideOfTriangle(*triangle, eq, p1)) {
      // We are modifying the constraint maybe it would be better to
      // not change the given constraint and just keep a variable for the new constraint
      edge_event_->constrained_edge.q = p1;
      triangle = triangle->NeighborAcross(point);
      EdgeEvent(ep, p1, triangle, p1);
    } else {
      HandleError("EdgeEvent - collinear points not supported");
    }
    return;
  }

  const Point* p2 = triangle->PointCW(point);
  const Orientation o2 = Orient2d(*eq, *p2, *ep);
  TRACE_OUT << "EdgeEvent - o1=" << o1 << "; o2=" << o2 << std::endl;
  if (o2 == COLLINEAR) {
    if (SetConstrainedEdgeIfSideOfTriangle(*triangle, eq, p2)) {
      // We are modifying the constraint maybe it would be better to
      // not change the given constraint and just keep a variable for the new constraint
      edge_event_->constrained_edge.q = p2;
      triangle = triangle->NeighborAcross(point);
      EdgeEvent(ep, p2, triangle, p2);
    } else {
      HandleError("EdgeEvent - collinear points not supported");
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

Node& Sweep::NewFrontTriangle(const Point* point, Node& node)
{
  Triangle* triangle = tcx_.AddTriangle(point, node.point, node.next->point);

  TRACE_OUT << "NewFrontTriangle - triangle=" << *triangle << std::endl;

  int i, oi;
  assert(node.triangle);
  triangle->MarkNeighbor(*node.triangle, i, oi);
  triangle->SetConstrainedEdge(i, node.triangle->IsConstrainedEdge(oi));

  Node* new_node = NewNode(point);

  new_node->next = node.next;
  new_node->prev = &node;
  node.next->prev = new_node;
  node.next = new_node;

  node.ResetTriangle();
  node.SetTriangle(triangle);
  new_node->SetTriangle(triangle);

  LegalizePush(*triangle);

  return *new_node;
}

void Sweep::Fill(Node** node)
{
  assert(node); assert(*node);
  Node* filled_node = *node;
  Triangle* triangle = tcx_.AddTriangle(filled_node->prev->point, filled_node->point, filled_node->next->point);

  TRACE_OUT << "Fill - triangle=" << *triangle << std::endl;

  if (filled_node->prev->triangle) {
    int i, oi;
    triangle->MarkNeighbor(*filled_node->prev->triangle, i, oi);
    triangle->SetConstrainedEdge(i, filled_node->prev->triangle->IsConstrainedEdge(oi));
  }
  if (filled_node->triangle) {
    int i, oi;
    triangle->MarkNeighbor(*filled_node->triangle, i, oi);
    triangle->SetConstrainedEdge(i, filled_node->triangle->IsConstrainedEdge(oi));
  }

  // Update the front
  assert(front_);
  Node* prev_node = front_->RemoveNode(*node);
  FreeNode(*node);
  prev_node->ResetTriangle();
  prev_node->SetTriangle(triangle);
  *node = prev_node;

  LegalizePush(*triangle);
}

namespace {

// Compute the angle formed by half lines OA and OB
Angle OABAngle(const Point* origin, const Point* pa, const Point* pb)
{
  /* In the complex plane: a = ax + i*ay; b = bx + i*by
   *
   * The angle is the argument of the complex number b/a:
   * b/a = (bx + by*i)/(ax + ay*i) = [(ax*bx + ay*by) + i*(ax*by - ay*bx)] / norm(a)^2
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
  const double x = ax * bx + ay * by;
  const double y = ax * by - ay * bx;
  return Angle(x, y);
}

// In [Zalik 2005] the criteria to fill a hole after a point event is to find
// if the angle the node forms with its neighbors is less that PI/2
// This can be evaluated like this:
//   const Angle node_angle = OABAngle(node->point, node->next->point, node->prev->point);
//   return node_angle.quadrant() == Angle::Quadrant::ONE;
//
// In the following function the criteria is refined. In case the hole is fillable but does
// not meet the criteria above, we also test the neighbor nodes on the same criteria,
// assuming the initial node was filled.
bool NodeFillCriteria(const Node* node, double fill_direction)
{
  // Copy of the NodeFillCriteria
  const Angle node_angle = OABAngle(node->point, node->next->point, node->prev->point);

  switch (node_angle.quadrant())
  {
    case Angle::Quadrant::ONE:
      return true;

    case Angle::Quadrant::TWO:
      // The hole could be filled. Assuming it is filled, evaluate if it the next neighbor node
      // (in the fill direction) will meet the fill hole criteria
      if (std::signbit(fill_direction) && node->prev->prev) {
        const Angle left_fill_angle = OABAngle(node->prev->point, node->next->point, node->prev->prev->point);
        return left_fill_angle.quadrant() == Angle::Quadrant::ONE;
      }
      if (!std::signbit(fill_direction) && node->next->next) {
        const Angle right_fill_angle = OABAngle(node->next->point, node->next->next->point, node->prev->point);
        return right_fill_angle.quadrant() == Angle::Quadrant::ONE;
      }
      return false;

    case Angle::Quadrant::THREE:
      // Not fillable
      return false;

    case Angle::Quadrant::FOUR:
    default:
      assert(0);
      return false;
  }
}

// The basin angle is decided against the horizontal line
bool LeftBasinCriteria(const Node& node)
{
  static const double TAN_PI_DIV_4 = std::tan(M_PI / 4.0);
  assert(node.prev && node.prev->prev);
  const Point v = *node.point - *node.prev->prev->point;
  const Angle slope(v.x, v.y);
  return slope.quadrant() == Angle::Quadrant::ONE && slope.tan() > TAN_PI_DIV_4;
}

// The basin angle is decided against the horizontal line
bool RightBasinCriteria(const Node& node)
{
  static const double TAN_3_PI_DIV_4 = std::tan(3.0 * M_PI / 4.0);
  assert(node.next && node.next->next);
  const Point v = *node.point - *node.next->next->point;
  const Angle slope(v.x, v.y);
  return slope.quadrant() == Angle::Quadrant::TWO && slope.tan() < TAN_3_PI_DIV_4;
}

} // namespace

void Sweep::FillAdvancingFront(Node& n)
{
  // Fill right holes
  Node* node = n.next;
  while (node && node->next && NodeFillCriteria(node, RIGHT_DIRECTION)) {
    Fill(&node);
    node = node->next;
  }

  // Fill left holes
  node = n.prev;
  while (node && node->prev && NodeFillCriteria(node, LEFT_DIRECTION)) {
    Fill(&node);
    // node is set to node->prev in Fill
  }

  // Fill right basin
  if (n.next && n.next->next && RightBasinCriteria(n)) {
    FillRightBasin(n);
  }

  // Fill left basin
  if (n.prev && n.prev->prev && LeftBasinCriteria(n)) {
    FillLeftBasin(n);
  }
}

void Sweep::LegalizePush(Triangle& triangle)
{
  legalize_stack_.emplace(&triangle, 0u);
}

void Sweep::Legalize() {
  std::size_t nb_flips = 0;
  const std::size_t max_nb_flips = 4 * tcx_.points_.size() * tcx_.points_.size();
  while (!legalize_stack_.empty()) {
    if (nb_flips > max_nb_flips) {
      // Note that this criteria for detecting an infinite loop is completely arbitrary!
      // We set a maximum number of flips equal to the square of the maximum number of triangles for this point set.
      // Actually the maximum number of triangulations (i.e. the space explored by the Legalization process) is
      // bounded by the Catalan numbers, which grows exponentially.
      HandleError("Legalize - Potential infinite loop was catched in Legalize");
    }
    const PendingLegalization legalize = legalize_stack_.top();
    legalize_stack_.pop();
    info_.max_legalize_depth = std::max(info_.max_legalize_depth, legalize.depth);
    Triangle* t = legalize.triangle;

    // To legalize a triangle we start by finding if any of the three edges violate the Delaunay condition
    for (int i = 0; i < 3; i++) {
      if (t->IsDelaunayEdge(i))
        continue;
      if (t->IsConstrainedEdge(i))
        continue;

      Triangle* ot = t->GetNeighbor(i);

      if (ot) {
        const Point* p = t->GetPoint(i);
        const Point* op = ot->OppositePoint(*t, p);
        int oi = ot->Index(op);
        assert(0 <= oi && oi < 3);
        assert(!ot->IsConstrainedEdge(oi));

        bool inside = InCircle(*p, *t->GetPoint((i + 1) % 3), *t->GetPoint((i + 2) % 3), *op);
        if (inside) {
          // If the Delaunay criteria is not verified after the triangle flip
          // the legalization process could go into an infinite loop.
          // This would always ne a consequence of numerical precision issues in the geometric predicates.
          assert(!InCircle(*t->GetPoint((i + 2) % 3), *p, *op, *t->GetPoint((i + 1) % 3)));

          // Lets rotate shared edge one vertex CW to legalize it (create a Delaunay pair)
          RotateTrianglePair(*t, i, *ot, oi, true);
          nb_flips++;
          TRACE_OUT << "Legalized - depth=" << legalize.depth << "; t=" << *t << "; ot=" << *ot << std::endl;

          // We now got one valid Delaunay Edge shared by two triangles
          // This gives us 4 new edges to check for Delaunay: This function is called recursively
          legalize_stack_.emplace(t, legalize.depth + 1);
          legalize_stack_.emplace(ot, legalize.depth + 1);

          // If triangle have been legalized no need to check the other edges since
          // the recursive legalization will handle those so we can end here
          break;
        } else {
          // The shared edge is Delaunay
          t->SetDelaunayEdge(i, true);
          ot->SetDelaunayEdge(oi, true);
        }
      }
    }
  }
  info_.nb_triangle_flips += static_cast<unsigned int>(nb_flips);
}

struct Basin {
  Node* left_node;
  Node* bottom_node;
  Node* right_node;

  static constexpr double shallow_height_to_width_ratio = 0.5;

  Basin() :
    left_node(nullptr),
    bottom_node(nullptr),
    right_node(nullptr)
  {
  }

  bool IsShallow() const
  {
    assert(left_node); assert(bottom_node); assert(right_node);
    const double width = right_node->point->x - left_node->point->x;
    const double height = std::max(left_node->point->y, right_node->point->y) - bottom_node->point->y;
    return shallow_height_to_width_ratio * width > height;
  }
};

void Sweep::FillRightBasin(Node& node)
{
  Basin basin;

  if (Orient2d(*node.point, *node.next->point, *node.next->next->point) == CCW) {
    basin.left_node = node.next->next;
  } else {
    basin.left_node = node.next;
  }

  // Find the bottom node
  basin.bottom_node = basin.left_node;
  while (basin.bottom_node->next
         && basin.bottom_node->point->y >= basin.bottom_node->next->point->y) {
    basin.bottom_node = basin.bottom_node->next;
  }
  if (basin.bottom_node == basin.left_node) {
    // No valid basin
    return;
  }

  // Find the right node
  basin.right_node = basin.bottom_node;
  while (basin.right_node->next
         && basin.right_node->point->y < basin.right_node->next->point->y) {
    basin.right_node = basin.right_node->next;
  }
  if (basin.right_node == basin.bottom_node) {
    // No valid basin
    return;
  }

  FillBasin(basin);
}

void Sweep::FillLeftBasin(Node& node)
{
  Basin basin;

  if (Orient2d(*node.point, *node.prev->point, *node.prev->prev->point) == CW) {
    basin.right_node = node.prev->prev;
  } else {
    basin.right_node = node.prev;
  }

  // Find the bottom node
  basin.bottom_node = basin.right_node;
  while (basin.bottom_node->prev
         && basin.bottom_node->point->y >= basin.bottom_node->prev->point->y) {
    basin.bottom_node = basin.bottom_node->prev;
  }
  if (basin.bottom_node == basin.right_node) {
    // No valid basin
    return;
  }

  // Find the left node
  basin.left_node = basin.bottom_node;
  while (basin.left_node->prev
         && basin.left_node->point->y < basin.left_node->prev->point->y) {
    basin.left_node = basin.left_node->prev;
  }
  if (basin.left_node == basin.bottom_node) {
    // No valid basin
    return;
  }

  FillBasin(basin);
}

void Sweep::FillBasin(Basin& basin)
{
  if (basin.IsShallow())
    return;

  TRACE_OUT << "FillBasin - left_node->point=" << *basin.left_node->point << "; bottom_node->point=" << *basin.bottom_node->point << "; right_node->point=" << *basin.right_node->point << std::endl;

  Node*& bottom = basin.bottom_node;
  assert(bottom);
  while (bottom != basin.left_node && bottom != basin.right_node) {
    Fill(&bottom);
    bottom = (bottom->point->y >= bottom->next->point->y) ? bottom->next : bottom;
  }
}

void Sweep::FillEdgeEvent(const Edge* edge, Node* node)
{
  assert(edge_event_);
  if (edge_event_->right) {
    FillRightAboveEdgeEvent(edge, node);
  } else {
    FillLeftAboveEdgeEvent(edge, node);
  }
}

void Sweep::FillRightAboveEdgeEvent(const Edge* edge, Node* node)
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

void Sweep::FillRightBelowEdgeEvent(const Edge* edge, Node& node)
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

void Sweep::FillRightConcaveEdgeEvent(const Edge* edge, Node& node)
{
  {
    Node* next = node.next;
    assert(next);
    Fill(&next);
  }
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

void Sweep::FillRightConvexEdgeEvent(const Edge* edge, Node& node)
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

void Sweep::FillLeftAboveEdgeEvent(const Edge* edge, Node* node)
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

void Sweep::FillLeftBelowEdgeEvent(const Edge* edge, Node& node)
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

void Sweep::FillLeftConvexEdgeEvent(const Edge* edge, Node& node)
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

void Sweep::FillLeftConcaveEdgeEvent(const Edge* edge, Node& node)
{
  {
    Node* prev = node.prev;
    assert(prev);
    Fill(&prev);
  }
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
  assert(edge_event_);
  TRACE_OUT << "FlipEdgeEvent - "
            << "edge={ ep=" << *ep << ", eq=" << *eq << " }; "
            << "triangle=" << *t << "; "
            << "p=" << *p
            << std::endl;

  Triangle* ot = t->NeighborAcross(p);
  if (ot == nullptr)
  {
    HandleError("FlipEdgeEvent - null neighbor across");
  }
  assert(ot);
  const Point* op = ot->OppositePoint(*t, p);

  const auto index_p = t->Index(p);
  const auto index_op = ot->Index(op);

  if (InScanArea(*p, *t->GetPoint((index_p + 1) % 3),*t->GetPoint((index_p + 2) % 3), *op)) {
    // Lets rotate shared edge one vertex CW
    RotateTrianglePair(*t, index_p, *ot, index_op);
    info_.nb_triangle_flips++;

    if (p == eq && op == ep) {
      if (eq == edge_event_->constrained_edge.q && ep == edge_event_->constrained_edge.p) {
        t->SetConstrainedEdge(ep, eq);
        ot->SetConstrainedEdge(ep, eq);
        LegalizePush(*t);
        LegalizePush(*ot);
      } else {
        // TODO: I think one of the triangles should be legalized here?
      }
    } else {
      // Select next flip triangle
      if (Orient2d(*eq, *op, *ep) == CW) {
        std::swap(t, ot);
      }
      // t is the next flip triangle
      // ot is not crossing the edge {ep, eq} and will be legalized in post-order
      FlipEdgeEvent(ep, eq, t, p);
      LegalizePush(*ot);
    }
  } else {
    const Point* new_p = NextFlipPoint(ep, eq, *ot, op);
    FlipScanEdgeEvent(ep, eq, *t, *ot, new_p);
    EdgeEvent(ep, eq, t, p);
  }
}

const Point* Sweep::NextFlipPoint(const Point* ep, const Point* eq, Triangle& ot, const Point* op)
{
  const Point* result = nullptr;
  const Orientation o2d = Orient2d(*eq, *op, *ep);
  if (o2d == CW) {
    // Right
    result = ot.PointCCW(op);
  } else if (o2d == CCW) {
    // Left
    result = ot.PointCW(op);
  } else {
    assert(o2d == COLLINEAR);
    HandleError("[Unsupported] Opposing point on constrained edge");
  }
  return result;
}

void Sweep::FlipScanEdgeEvent(const Point* ep, const Point* eq, Triangle& flip_triangle,
                              Triangle& t, const Point* p)
{
  TRACE_OUT << "FlipScanEdgeEvent - "
            << "edge={ ep=" << *ep << ", eq=" << *eq << " }; "
            << "triangle=" << t << "; "
            << "p=" << *p << "; "
            << "flip_triangle=" << flip_triangle
            << std::endl;
  Triangle* ot_ptr = t.NeighborAcross(p);
  if (ot_ptr == nullptr) {
    HandleError("FlipScanEdgeEvent - null neighbor across");
  }
  assert(ot_ptr);

  const Point* op = ot_ptr->OppositePoint(t, p);
  if (op == nullptr) {
    HandleError("FlipScanEdgeEvent - null opposing point");
  }
  assert(op);

  const Point* p1 = flip_triangle.PointCCW(eq);
  const Point* p2 = flip_triangle.PointCW(eq);
  if (p1 == nullptr || p2 == nullptr) {
    HandleError("FlipScanEdgeEvent - null on either of points");
  }
  assert(p1); assert(p2);

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

void Sweep::HandleError(std::string_view msg)
{
  DeleteFront();
  while (!legalize_stack_.empty()) { legalize_stack_.pop(); }
  throw std::runtime_error(msg.data());
}

Node* Sweep::NewNode(const Point* p, Triangle* t)
{
  assert(p);
  Node* new_node = node_storage_->NewNode();
  *new_node = Node(p, t);
  if (t) { t->SetNode(*new_node); }
  return new_node;
}

void Sweep::FreeNode(Node* node)
{
  assert(node);
  node_storage_->DiscardNode(node);
}

} // namespace p2t
