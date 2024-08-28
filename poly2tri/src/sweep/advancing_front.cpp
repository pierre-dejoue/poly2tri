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

Node* AdvancingFront::RemoveNode(Node* deleted_node)
{
  assert(deleted_node);
  assert(deleted_node);
  assert(deleted_node->prev);   // The node is NOT the head of the list
  assert(deleted_node->next);   // The node is NOT the tail of the list
  deleted_node->prev->next = deleted_node->next;
  deleted_node->next->prev = deleted_node->prev;
  Node* prev_node = deleted_node->prev;
  deleted_node->prev = nullptr;
  deleted_node->ResetTriangle();
  if (search_node_ == deleted_node) { search_node_ = prev_node; }
  assert(prev_node);
  return prev_node;
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
