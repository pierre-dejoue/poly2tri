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

#include <poly2tri/common/shapes.h>

#include <cassert>
#include <ostream>
#include <utility>


namespace p2t {

// Advancing front node
struct Node {

  const Point* point;
  Triangle* triangle;     // This triangle (if set) should contain the edge between this node and the next one.

  Node* next;
  Node* prev;

  double value;

  Node(const Point* p, Triangle* t = nullptr) :
    point(p),
    triangle(t),
    next(nullptr),
    prev(nullptr),
    value()
  {
    assert(p);
    value = p->x;
  }

};

std::ostream& operator<<(std::ostream& out, const Node& node);

// Advancing front
class AdvancingFront {
public:

  AdvancingFront(Node& head, Node& tail);
  ~AdvancingFront();

  Node* head() const { return head_; }
  Node* tail() const { return tail_; }

  /// Locate insertion point along advancing front
  Node* LocateNode(double x);

  Node* LocatePoint(const Point* point);

private:

  Node* const head_;
  Node* const tail_;
  Node* search_node_;

  Node* FindSearchNode(double x);

};

std::ostream& operator<<(std::ostream& out, const AdvancingFront& front);

// Return a range (begin, end) of the inner nodes of the advancing front
// That is, the range of all nodes from head()->next->next to tail()->prev->prev
std::pair<Node*, Node*> GetInnerRange(AdvancingFront& front);

}
