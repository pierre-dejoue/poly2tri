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

#include "advancing_front.h"

#include "../common/utils.h"

#include <cassert>

namespace p2t {

AdvancingFront::AdvancingFront(Node& head, Node& tail) :
  head_(&head),
  tail_(&tail),
  search_node_(&head)
{
}

std::ostream& operator<<(std::ostream& out, const Node& node)
{
  out << "{ point=" << *node.point << "; triangle=";
  if (node.triangle) {
    out << *node.triangle;
    if (!node.triangle->IsInterior()) {
       out << "; EXTERIOR";
    }
  } else {
    out << "NULL";
  }
  return out << " }";
}

AdvancingFront::~AdvancingFront() = default;

Node* AdvancingFront::LocateNode(double x)
{
  Node* node = search_node_;

  if (x < node->value) {
    while ((node = node->prev) != nullptr) {
      if (x >= node->value) {
        search_node_ = node;
        return node;
      }
    }
  } else {
    while ((node = node->next) != nullptr) {
      if (x < node->value) {
        search_node_ = node->prev;
        return node->prev;
      }
    }
  }
  return nullptr;
}

Node* AdvancingFront::FindSearchNode(double x)
{
  UNUSED(x);

  // TODO: implement BST index
  return search_node_;
}

Node* AdvancingFront::LocatePoint(const Point* point)
{
  const double x_dir = (head_->point->x < tail_->point->x) ? 1.0 : -1.0;
  const double px = point->x;
  Node* node = FindSearchNode(px);
  const double nx = node->point->x;

  if (px == nx) {
    if (point != node->point) {
      // We might have two nodes with same x value for a short time
      if (point == node->prev->point) {
        node = node->prev;
      } else if (point == node->next->point) {
        node = node->next;
      } else {
        assert(0);
      }
    }
  } else if (x_dir * px < x_dir * nx) {
    while ((node = node->prev) != nullptr) {
      if (point == node->point)
        break;
    }
  } else {
    while ((node = node->next) != nullptr) {
      if (point == node->point)
        break;
    }
  }
  if(node) { search_node_ = node; }
  return node;
}

void AdvancingFront::MapTriangleToNodes(Triangle& t)
{
  for (int i = 0; i < 3; i++) {
    if (!t.GetNeighbor(i)) {
      Node* node = LocatePoint(t.PointCW(t.GetPoint(i)));
      if (node) { node->triangle = &t; }
    }
  }
}

std::ostream& operator<<(std::ostream& out, const AdvancingFront& front)
{
  const Node* head = front.head();
  const Node* tail = front.tail();
  assert(head != nullptr); assert(tail != nullptr);
  out << "{\n    head=" << *head << " -->\n";
  const Node* node = head->next;
  while (node != tail) {
    out << "         " << *node << " -->\n";
    node = node->next;
  }
  return out  << "    tail=" << *tail << " }";
}

std::pair<Node*, Node*> GetInnerRange(AdvancingFront& front)
{
  std::size_t idx = 1;
  Node* node = front.head()->next;
  Node* const tail = front.tail();
  Node* begin_node = nullptr;
  while (node != nullptr && node != tail && idx < 3)
  {
    idx++;
    node = node->next;
    if (idx == 2)
      begin_node = node;
  }
  if (idx < 3)
  {
    // idx = 1: h -> t
    // idx = 2: h -> n -> t
    return std::make_pair(front.head(), front.head());     // Return an empty range
  }
  // Else, the advancing range has a length at least equal to 4 (The shortest case is: h -> n1 -> n2 -> t, with begin_node = n2)
  assert(begin_node != nullptr);
  assert(tail->prev != nullptr);
  return std::make_pair(begin_node, tail->prev);
}

} // namespace p2t
