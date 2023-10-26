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

#include <poly2tri/common/orientation.h>

#include <cassert>
#include <cmath>

/**
 * Definition of the geometric predicates used in poly2tri
 *
 * The input is a generic Point type which can be as simple as:
 *
 *   struct P {
 *    using scalar = double;
 *    scalar x, y;
 *   };
 */
namespace p2t {

/**
 * <b>Determines the orientation of triangle ABC</b><br>
 * Uses the formula to calculate the signed area:<br>
 * <pre>
 * A[P1,P2,P3]  =  (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + (x3*y1 - y3*x1)
 *              =  (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3)
 * </pre>
 * The area is:
 *   Positive if the triangle is CCW
 *   Negative if it is CW
 *   Zero if it is degenerated.
 */
template <typename P>
inline Orientation Orient2d(const P& pa, const P& pb, const P& pc)
{
  using F = typename P::scalar;

  F oabc = (pa.x - pc.x) * (pb.y - pc.y) - (pa.y - pc.y) * (pb.x - pc.x);
  // Using a tolerance here fails on concave-by-subepsilon boundaries
  // e.g.  if (oabc > -EPSILON && oabc < EPSILON) {
  // Using == on double makes -Wfloat-equal warnings yell at us
  if (std::fpclassify(oabc) == FP_ZERO) {
    return COLLINEAR;
  } else if (oabc > 0) {
    return CCW;
  }
  return CW;
}


/**
 * <b>Determines if d is inside the circumcircle of triangle abc</b><br>
 * <b>Requirement</b>:<br>
 * 1. a,b and c form a triangle.<br>
 * 2. a and d are known to be on opposite side of bc<br>
 * <pre>
 *                a
 *                +
 *               / \
 *              /   \
 *            b/     \c
 *            +-------+
 *           /    d    \
 *          /     +     \
 * </pre>
 * <b>Fact</b>: d has to be in area A delimited by semi-lines ab and ac to have a chance to be inside
 * the circumcircle of triangle abc<br>
 *  d is outside A if orient2d(a,b,d) or orient2d(c,a,d) is CW<br>
 *  This preknowledge gives us a way to optimize the incircle test
 *
 * @param pa - triangle point, opposite to d
 * @param pb - triangle point
 * @param pc - triangle point
 * @param pd - point opposite to a
 * @return true if d is inside the circumcircle of triangle abc, false if d is on or outside that circle
 */
template <typename P>
inline bool InCircle(const P& pa, const P& pb, const P& pc, const P& pd)
{
  using F = typename P::scalar;

  const F adx = pa.x - pd.x;
  const F ady = pa.y - pd.y;
  const F bdx = pb.x - pd.x;
  const F bdy = pb.y - pd.y;

  const F adxbdy = adx * bdy;
  const F bdxady = bdx * ady;
  const F oabd = adxbdy - bdxady;

  if (oabd <= 0)
    return false;

  const F cdx = pc.x - pd.x;
  const F cdy = pc.y - pd.y;

  const F cdxady = cdx * ady;
  const F adxcdy = adx * cdy;
  const F ocad = cdxady - adxcdy;

  if (ocad <= 0)
    return false;

  const F bdxcdy = bdx * cdy;
  const F cdxbdy = cdx * bdy;
  const F obcd = bdxcdy - cdxbdy;

  const F adsqnorm = adx * adx + ady * ady;
  const F bdsqnorm = bdx * bdx + bdy * bdy;
  const F cdsqnorm = cdx * cdx + cdy * cdy;

  const F det = adsqnorm * obcd + bdsqnorm * ocad + cdsqnorm * oabd;

  return det > 0;
}

} // namespace p2t
