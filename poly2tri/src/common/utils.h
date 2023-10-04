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

// Otherwise #defines like M_PI are undeclared under Visual Studio
#define _USE_MATH_DEFINES

#include <poly2tri/common/shapes.h>

#include <cmath>
#include <limits>

// C99 removes M_PI from math.h
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

// For unused function arguments
#define UNUSED(x) (void)(x)

namespace p2t {

constexpr double EPSILON = 1e-12;

constexpr double LEFT_DIRECTION = -1.0;
constexpr double RIGHT_DIRECTION = 1.0;

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
inline Orientation Orient2d(const Point& pa, const Point& pb, const Point& pc)
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
inline bool InScanArea(const Point& pa, const Point& pb, const Point& pc, const Point& pd)
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

class Angle {
public:
  // Quadrants: ONE = [0, PI/2), TWO = [PI/2, PI), THREE = [PI, 3*PI/2), FOUR = [3*PI/2, 2*PI)
  enum class Quadrant {
    ONE = 0,
    TWO = 1,
    FOUR = 2,
    THREE = 3
  };

  Angle(double x, double y) : x_(x), y_(y), q_(Quadrant::ONE)
  {
    const unsigned int q = (std::signbit(x) ? 1u : 0u) ^ (std::signbit(y) ? 2u : 0u);
    q_ = static_cast<Quadrant>(q);
  }

  Quadrant quadrant() const
  {
    return q_;
  }

  double tan() const
  {
    if (std::fpclassify(x_) == FP_ZERO) {
      return (std::signbit(y_) ? -1.0 : 1.0) * std::numeric_limits<double>::max();
    } else {
      return y_ / x_;
    }
  }

private:
  double x_, y_;
  Quadrant q_;
};

}
