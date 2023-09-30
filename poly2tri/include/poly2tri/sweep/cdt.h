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

#include "../common/dll_symbol.h"
#include "../common/shapes.h"
#include "policy.h"

#include <memory>
#include <vector>

/**
 *
 * @author Mason Green <mason.green@gmail.com>
 *
 */

namespace p2t {

struct Point;
struct Edge;
class Sweep;
class SweepContext;
class Triangle;

class P2T_DLL_SYMBOL CDT
{
public:

  /**
   * Default constructor
   */
  CDT();

  /**
   * Constructor - add polyline with non repeating points
   *
   * Note: prefer using the default constructor followed by a call to AddPolyline()
   *
   * @param polyline
   */
  CDT(const std::vector<Point*>& polyline);   // Kept for backward compatibility

  /**
   * Destructor - clean up memory
   */
  ~CDT();

  /**
   * Add the outer polyline
   *
   *  - Must be unique
   *  - Must have no repeating points
   *  - It is closed (the first and last points form a constrained edge)
   *  - Call this before the other methods: AddHole, AddPoint
   *
   * @param polyline
   * @param num_points
   * @param stride (in bytes) If zero, assume a contiguous array of Points
   */
  void AddPolyline(const std::vector<Point*>& polyline);  // Kept for backward compatibility
  void AddPolyline(const Point* const* polyline, std::size_t num_points);
  void AddPolyline(const Point* polyline, std::size_t num_points, std::size_t stride = 0u);

  /**
   * Add a hole
   *
   * @param polyline
   * @param num_points
   * @param stride (in bytes) If zero, assume a contiguous array of Points
   */
  void AddHole(const std::vector<Point*>& polyline);      // Kept for backward compatibility
  void AddHole(const Point* const* polyline, std::size_t num_points);
  void AddHole(const Point* polyline, std::size_t num_points, std::size_t stride = 0u);

  /**
   * Add a steiner point
   *
   * @param point
   */
  void AddPoint(const Point* point);

  /**
   * Add several steiner points
   *
   * @param points
   * @param num_points
   * @param stride (in bytes) If zero, assume a contiguous array of Points
   */
  void AddPoints(const std::vector<Point*>& points);      // Kept for backward compatibility
  void AddPoints(const Point* const* points, std::size_t num_points);
  void AddPoints(const Point* points, std::size_t num_points, std::size_t stride = 0u);

  /**
   * Triangulate - do this AFTER you've added the polyline, holes, and Steiner points
   */
  void Triangulate(Policy triangulation_policy = Policy::OuterPolygon);

  /**
   * Get the number of CDT triangles
   */
  std::size_t GetTrianglesCount() const;

  /**
   * Get CDT triangles
   *
   * OutIt: An output iterator of Triangle*
   */
  const std::vector<std::unique_ptr<Triangle>>& GetTriangles() const;
  template <typename OutIt>
  void GetTriangles(OutIt triangle_pointer_dest) const;

  /**
   * Get the number of points that were given as input by the user
   */
  std::size_t GetInputPointsCount() const;

  /**
   * Get the list of all points that were given as input by the user
   *
   * For test and debug purpose only.
   */
  std::vector<const Point*> GetInputPoints() const;

  /**
   * Get the number of edges that were given as input by the user
   */
  std::size_t GetInputEdgesCount() const;

  /**
   * Get the list of all edges that were given as input by the user
   *
   * For test and debug purpose only.
   */
  std::vector<Edge*> GetInputEdges() const;

private:

  std::unique_ptr<SweepContext> sweep_context_;

};

// Utility function to help the transition from the original API: std::vector<p2t::Triangle*> CDT::GetTriangles();
std::vector<Triangle*> GetTrianglesAsVector(const p2t::CDT& cdt);

//
// Implementations
//

template <typename OutIt>
void CDT::GetTriangles(OutIt triangle_pointer_dest) const
{
  for (const auto& t : GetTriangles())
  {
    *triangle_pointer_dest++ = t.get();
  }
}

}
