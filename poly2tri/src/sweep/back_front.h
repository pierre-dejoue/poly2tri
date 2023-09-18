/*
 * Poly2Tri Copyright (c) 2023, Poly2Tri Contributors
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

#include "advancing_front.h"

#include <cassert>
#include <cstddef>

namespace p2t
{

class Sweep;

using BackFront = AdvancingFront;

/**
 * Traverse all the back triangles, meaning the triangles which contains at least one of the two artificial points (head and tail.)
 * The traversal starts with the tail and ends with the head.
 *
 * The back front is made of three families of triangles, encountered in the following order when traversing from tail (T) to head (H):
 * - A fan of triangles with T as the pivot point
 * - A unique triangle which contains the segment HT (we call it the "middle" triangle)
 *-  A fan of triangles with H as the pivot point
 *
 * An action on a back triangle is a function with prototype void(Triangle* triangle, int pivot, bool mid, bool last), with argumentss:
 *  - triangle: pointer to the back triangle
 *  - pivot (= 0, 1 or 2): the index of the "pivot" vertex in the triangle, which is either the head or the tail
 *  - mid (boolean): true iff this is the "middle" triangle
 *  - last (boolean): true iff this is the last triangle of the traversal
 */
template <typename ACTION>
void TraverseBackTriangles(AdvancingFront& adv_front, ACTION action)
{
  const Node* head = adv_front.head();
  const Node* tail = adv_front.tail();
  const Node* prev_tail = tail->prev;
  assert(prev_tail != nullptr);
  Triangle* triangle = prev_tail->triangle;
  assert(triangle != nullptr);
  Triangle* neighbor = nullptr;
  const Point* pivot = tail->point;
  while ((neighbor = triangle->NeighborCCW(pivot)) != nullptr) { triangle = neighbor; }
  while (triangle != nullptr)
  {
    bool mid = false;
    int index = triangle->Index(pivot);
    Triangle* next = triangle->NeighborCW(pivot);
    if (next == nullptr && pivot == tail->point)
    {
      mid = true;
      assert(triangle->PointCW(pivot) == head->point);
      pivot = head->point;
      next = triangle->NeighborCW(pivot);
    }
    action(triangle, index, mid, next == nullptr);
    triangle = next;
  }
}

}
