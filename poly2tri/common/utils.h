/*
 * Poly2Tri Copyright (c) 2009-2018, Poly2Tri Contributors
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

#ifndef UTILS_H
#define UTILS_H

// Otherwise #defines like M_PI are undeclared under Visual Studio
#define _USE_MATH_DEFINES

#include "shapes.h"

#include <cmath>
#include <exception>

// C99 removes M_PI from math.h
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace p2t {

const double PI_3div4 = 3 * M_PI / 4;
const double PI_div2 = 1.57079632679489661923;
const double EPSILON = 1e-12;

enum Orientation { CW, CCW, COLLINEAR };

/**
 * <b>Determines the orientation of triangle abc</b><br>
 * Uses the formula to calculate the signed area:<br>
 * <pre>
 * A[P1,P2,P3]  =  (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + (x3*y1 - y3*x1)
 *              =  (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3)
 * </pre>
 * Positive if CCW<br>
 * Negative if CW<br>
 * 0 if collinear<br>
 */
Orientation Orient2d(const Point& pa, const Point& pb, const Point& pc)
{
  double oabc = (pa.x - pc.x) * (pb.y - pc.y) - (pa.y - pc.y) * (pb.x - pc.x);
// Using a tolerance here fails on concave-by-subepsilon boundaries
//   if (oabc > -EPSILON && oabc < EPSILON) {
// Using == on double makes -Wfloat-equal warnings yell at us
  if (std::fpclassify(oabc) == FP_ZERO) {
    return COLLINEAR;
  } else if (oabc > 0) {
    return CCW;
  }
  return CW;
}

/**
 * <b>Determines if semi-line ad intersects triangle abc</b><br>
 * Corner case: if d is collinear with a and b, or a and c, then return false
 */
bool InScanArea(const Point& pa, const Point& pb, const Point& pc, const Point& pd)
{
  double odab = (pd.x - pb.x) * (pa.y - pb.y) - (pa.x - pb.x) * (pd.y - pb.y);
  if (odab <= EPSILON) {
    return false;
  }

  double oadc = (pa.x - pc.x) * (pd.y - pc.y) - (pd.x - pc.x) * (pa.y - pc.y);
  if (oadc <= EPSILON) {
    return false;
  }
  return true;
}

}

#endif
