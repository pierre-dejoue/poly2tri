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

#include <cassert>
#include <cmath>
#include <cstddef>
#include <ostream>

namespace p2t {

struct Point {

  double x, y;

  Point() : x(0.0), y(0.0) {}
  Point(double x_, double y_) : x(x_), y(y_) {}

  /// Set this point to all zeros.
  void set_zero()
  {
    x = 0.0;
    y = 0.0;
  }

  /// Set this point to some specified coordinates.
  void set(double x_, double y_)
  {
    x = x_;
    y = y_;
  }

  /// Negate this point.
  Point operator -() const
  {
    Point v;
    v.set(-x, -y);
    return v;
  }

  /// Add a point to this point.
  void operator +=(const Point& v)
  {
    x += v.x;
    y += v.y;
  }

  /// Subtract a point from this point.
  void operator -=(const Point& v)
  {
    x -= v.x;
    y -= v.y;
  }

  /// Multiply this point by a scalar.
  void operator *=(double a)
  {
    x *= a;
    y *= a;
  }

  /// Get the length of this point (the norm).
  double Length() const
  {
    return std::sqrt(x * x + y * y);
  }

  /// Convert this point into a unit point. Returns the Length.
  double Normalize()
  {
    const double len = Length();
    x /= len;
    y /= len;
    return len;
  }

};

inline bool cmp(const Point* a, const Point* b)
{
  if (a->y < b->y) {
    return true;
  } else if (a->y == b->y) {
    // Make sure q is point with greater x value
    if (a->x < b->x) {
      return true;
    }
  }
  return false;
}

/// Add two points_ component-wise.
inline Point operator +(const Point& a, const Point& b)
{
  return Point(a.x + b.x, a.y + b.y);
}

/// Subtract two points_ component-wise.
inline Point operator -(const Point& a, const Point& b)
{
  return Point(a.x - b.x, a.y - b.y);
}

/// Multiply point by scalar
inline Point operator *(double s, const Point& a)
{
  return Point(s * a.x, s * a.y);
}

inline bool operator ==(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}

inline bool operator !=(const Point& a, const Point& b)
{
  return !(a.x == b.x) || !(a.y == b.y);
}

/// Peform the dot product on two vectors.
inline double Dot(const Point& a, const Point& b)
{
  return a.x * b.x + a.y * b.y;
}

/// Perform the cross product on two vectors. In 2D this produces a scalar.
inline double Cross(const Point& a, const Point& b)
{
  return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a point and a scalar. In 2D this produces
/// a point.
inline Point Cross(const Point& a, double s)
{
  return Point(s * a.y, -s * a.x);
}

/// Perform the cross product on a scalar and a point. In 2D this produces
/// a point.
inline Point Cross(double s, const Point& a)
{
  return Point(-s * a.y, s * a.x);
}

}
