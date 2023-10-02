Poly2tri
========

## About this fork of the library

This fork of [jhasse/poly2tri](https://github.com/jhasse/poly2tri) was created by Pierre Dejoue [:email:](mailto:pierre.dejoue@gmail.com).

In this repository I share a number of changes I made to the library for my own purpose. Some of thoses changes impact the API and are listed below.

There is no intent to provide any level of maintainance or support of this repository. It is open source and closed to contributions. That being stated, feel free to [submit an issue](https://github.com/pierre-dejoue/poly2tri/issues/new) for a bug report, or a question.

## Foreword of the original poly2tri library

Since there are no Input validation of the data given for triangulation you need
to think about this. Poly2Tri does not support repeat points within epsilon.

* If you have a cyclic function that generates random points make sure you don't
  add the same coordinate twice.
* If you are given input and aren't sure same point exist twice you need to
  check for this yourself.
* Only simple polygons are supported. You may add holes or interior Steiner points
* Interior holes must not touch other holes, nor touch the polyline boundary
* Use the library in this order:
  1. Initialize CDT with a simple polyline (this defines the constrained edges)
  2. Add holes if necessary (also simple polylines)
  3. Add Steiner points
  4. Triangulate

Make sure you understand the preceding notice before posting an issue. If you have
an issue not covered by the above, include your data-set with the problem.
The only easy day was yesterday; have a nice day. <Mason Green>

Poly2tri Installation Guide
===========================

API Changes
-----------

Compared to the forked [repository](https://github.com/jhasse/poly2tri).

### Breaking

- Renamed method Triangle::CircumcircleContains to fix a typo in the name
- Renamed methods Triangle::{Get/Set}DelaunayEdge{CW/CCW} to fix a typo in their name and also change the 'Get' prefix to 'Is'.
- Removed methods Triangle::MarkConstrainedEdge and use Triangle::SetContrainedEdge instead. Those are very unlikely to be used in client code.
- Members Triangle::contrained_edge[] and Triangle::delaunay_edge[] are now private
- Removed function IsDelaunay
- Method Triangle::GetPoint(index) now returns a const Point*
- Many other methods in class Triangle have changed to take or return const Point* instead of references.
  Most of those methods, despite being public, are unlikely to be used in client code.
- Remove CDT::GetMap(). That method was for debug purpose only and imposed an extra memory cost.
- Method CDT::GetTriangles() returns a vector of unique_ptr by const ref. This is likely to break some client
  code. A free function p2t::GetTrianglesAsVector(cdt) has been added to facilitate the transition.
- Removed p2t::Triangle::DebugPrint(). Replaced by overloaded operator<<.

### Non-Breaking

- Orientation enum is now public
- Added method Triangle::GetOrientation()
- Added a default constructor to class CDT
- Added methods CDT::AddPolyline(polyline) and CDT::AddPoints(points)
- Make the Point structure safer and more flexible:
    - Remove the edge_list from Point
    - New API CDT::Add{Polyline/Hole/Points} with pointer + size arguments
    - Never modify the user's Point data (const Point*)
- Add CDT::GetTrianglesCount()
- Add CDT::GetTriangles(output_iterator)
- Add CDT::AddOpenPolyline(...)
- Neighbor triangles in the triangulation result are now consistent (no links to triangles exterior to the CDT)
- Add CDT::LastTriangulationInfo()

Dependencies
------------

Core poly2tri lib:

* C++17
* Standard Template Library (STL)

Unit tests:

* Boost (filesystem, test framework)

Testbed:

* OpenGL
* [GLFW](http://glfw.sf.net)

Build the library
-----------------

With the ninja build system installed:

```
mkdir build && cd build
cmake -GNinja ..
cmake --build .
```

Build and run with unit tests
----------------------------

With the ninja build system:

```
mkdir build && cd build
cmake -GNinja -DP2T_BUILD_TESTS=ON ..
cmake --build .
ctest --output-on-failure
```

Build with the testbed
----------------------

```
mkdir build && cd build
cmake -GNinja -DP2T_BUILD_TESTBED=ON ..
cmake --build .
```

Build with meson
----------------

```
meson setup build && cd build
meson compile
```

Running the Examples
--------------------

Load data points from a file:
```
build/testbed/p2t <filename> <center_x> <center_y> <zoom>
```
Load data points from a file and automatically fit the geometry to the window:
```
build/testbed/p2t <filename>
```
Random distribution of points inside a constrained box:
```
build/testbed/p2t random <num_points> <box_radius> <zoom>
```
Examples:
```
build/testbed/p2t testbed/data/dude.dat 350 500 3.0

build/testbed/p2t testbed/data/nazca_monkey.dat

build/testbed/p2t random 10   100 5.0
build/testbed/p2t random 1000 100 5.0
```

Testbed Controls
----------------

- Use the UP and DOWN arrow keys to zoom in/out
- Use BACKSPACE to reset zoom
- Use ESCAPE to quit

References
==========

- Domiter V. and Zalik B. (2008) Sweep‐line algorithm for constrained Delaunay triangulation
- Zalik B. (2005) An efficient sweep-line Delaunay triangulation algorithm
- https://gemma.feri.um.si/delaunay-triangulation/
- FlipScan by library author Thomas Åhlén:

![FlipScan](doc/FlipScan.png)
