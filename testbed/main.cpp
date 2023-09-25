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
#include <poly2tri/poly2tri.h>

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

bool ParseFile(std::string filename, std::vector<p2t::Point>& out_polyline,
               std::vector<std::vector<p2t::Point>>& out_holes, std::vector<p2t::Point>& out_steiner);
std::pair<p2t::Point, p2t::Point> BoundingBox(const std::vector<p2t::Point>& polyline);
void GenerateRandomPointDistribution(size_t num_points, double min, double max,
                                     std::vector<p2t::Point>& out_polyline,
                                     std::vector<std::vector<p2t::Point>>& out_holes,
                                     std::vector<p2t::Point>& out_steiner);
void Init(int window_width, int window_height);
void ShutDown(int return_code);
void MainLoop(double initial_zoom);
void Draw(const double zoom);
void ConstrainedColor(bool constrained);
double StringToDouble(const std::string& s);
double Random(double (*fun)(double), double xmin, double xmax);
double Fun(double x);

/// Default window size
constexpr int default_window_width = 800;
constexpr int default_window_height = 600;

/// Autozoom border (percentage)
const double autozoom_border = 0.05;

/// Flip Y axis
constexpr bool flag_flip_y = false;

/// Convex hull triangulation
constexpr bool convex_hull_triangulation = false;

/// Create a random distribution of points?
bool random_distribution = false;

/// Screen center
double cx = 0.0;
double cy = 0.0;

/// Constrained Delaunay triangles
std::vector<p2t::Triangle*> triangles;

/// Polylines
std::vector<p2t::Point> polyline;
std::vector<std::vector<p2t::Point>> holes;
std::vector<p2t::Point> steiner;

/// GLFW Window
GLFWwindow* window = nullptr;

int window_width = default_window_width;
int window_height = default_window_height;

int main(int argc, char* argv[])
{
  std::string filename;
  size_t num_points = 0u;
  double max, min;
  double zoom;

  if (argc != 2 && argc != 5) {
    std::cout << "-== USAGE ==-" << std::endl;
    std::cout << "Load Data File: p2t <filename> <center_x> <center_y> <zoom>" << std::endl;
    std::cout << "  Example: build/testbed/p2t testbed/data/dude.dat 350 500 3" << std::endl;
    std::cout << "Load Data File with Auto-Zoom: p2t <filename>" << std::endl;
    std::cout << "  Example: build/testbed/p2t testbed/data/nazca_monkey.dat" << std::endl;
    std::cout << "Generate Random Polygon: p2t random <num_points> <box_radius> <zoom>" << std::endl;
    std::cout << "  Example: build/testbed/p2t random 100 1 500" << std::endl;
    return 1;
  }

  // If true, adjust the zoom settings to fit the input geometry to the window
  const bool autozoom = (argc == 2);

  if (!autozoom && std::string(argv[1]) == "random") {
    num_points = atoi(argv[2]);
    random_distribution = true;
    char* pEnd;
    max = strtod(argv[3], &pEnd);
    min = -max;
    cx = cy = 0.0;
    zoom = atof(argv[4]);
  } else {
    filename = std::string(argv[1]);
    if (!autozoom) {
      cx = atof(argv[2]);
      cy = atof(argv[3]);
      zoom = atof(argv[4]);
    }
  }

  if (random_distribution) {
    GenerateRandomPointDistribution(num_points, min, max, polyline, holes, steiner);
  } else {
    // Load pointset from file
    if (!ParseFile(filename, polyline, holes, steiner)) {
      return 2;
    }
  }

  if (autozoom) {
    assert(0.0 <= autozoom_border && autozoom_border < 1.0);
    const auto bbox = BoundingBox(polyline);
    p2t::Point center = bbox.first + bbox.second;
    center *= 0.5;
    cx = center.x;
    cy = center.y;
    p2t::Point sides = bbox.second - bbox.first;
    zoom = 2.0 * (1.0 - autozoom_border) * std::min((double)default_window_width / sides.x, (double)default_window_height / sides.y);
    std::cout << "center_x = " << cx << std::endl;
    std::cout << "center_y = " << cy << std::endl;
    std::cout << "zoom = " << zoom << std::endl;
  }

  Init(default_window_width, default_window_height);

  /*
   * Perform triangulation!
   */

  double init_time = glfwGetTime();

  /*
   * STEP 1: Create CDT and add primary polyline
   * NOTE: polyline must be a simple polygon. The polyline's points
   * constitute constrained edges. No repeat points!!!
   */
  p2t::CDT cdt;

  /*
   * STEP 2: Add holes or Steiner points
   */
  if (convex_hull_triangulation) {
    cdt.AddPoints(polyline.data(), polyline.size());
    for (const auto& hole : holes) {
      assert(!hole.empty());
      cdt.AddPoints(hole.data(), hole.size());
    }
    cdt.AddPoints(steiner.data(), steiner.size());
  } else {
    cdt.AddPolyline(polyline.data(), polyline.size());
    for (const auto& hole : holes) {
      assert(!hole.empty());
      cdt.AddHole(hole.data(), hole.size());
    }
    cdt.AddPoints(steiner.data(), steiner.size());
  }

  /*
   * STEP 3: Triangulate!
   */
  const p2t::Policy policy = convex_hull_triangulation ? p2t::Policy::ConvexHull : p2t::Policy::OuterPolygon;
  cdt.Triangulate(policy);

  double dt = glfwGetTime() - init_time;

  triangles.reserve(cdt.GetTrianglesCount());
  cdt.GetTriangles(std::back_inserter(triangles));
  const size_t points_in_holes =
      std::accumulate(holes.cbegin(), holes.cend(), size_t(0),
                      [](size_t cumul, const std::vector<p2t::Point>& hole) { return cumul + hole.size(); });

  std::cout << "Number of primary constrained edges = " << polyline.size() << std::endl;
  std::cout << "Number of holes = " << holes.size() << std::endl;
  std::cout << "Number of constrained edges in holes = " << points_in_holes << std::endl;
  std::cout << "Number of Steiner points = " << steiner.size() << std::endl;
  std::cout << "Total number of points = " << (polyline.size() + points_in_holes + steiner.size()) << std::endl;
  std::cout << "Number of triangles = " << triangles.size() << std::endl;
  std::cout << "Elapsed time (ms) = " << dt * 1000.0 << std::endl;

  MainLoop(zoom);

  ShutDown(0);
  return 0;
}

bool ParseFile(std::string filename, std::vector<p2t::Point>& out_polyline, std::vector<std::vector<p2t::Point>>& out_holes,
               std::vector<p2t::Point>& out_steiner)
{
  enum ParserState {
    Polyline,
    Hole,
    Steiner,
  };
  ParserState state = Polyline;
  std::vector<p2t::Point>* hole = nullptr;
  try {
    std::string line;
    std::ifstream myfile(filename);
    if (myfile.is_open()) {
      while (!myfile.eof()) {
        getline(myfile, line);
        if (line.empty()) {
          break;
        }
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(tokens));
        if (tokens.empty()) {
          break;
        } else if (tokens.size() == 1u) {
          const auto token = tokens[0];
          if (token == "HOLE") {
            state = Hole;
            out_holes.emplace_back();
            hole = &out_holes.back();
          } else if (token == "STEINER") {
            state = Steiner;
          } else {
            throw std::runtime_error("Invalid token [" + token + "]");
          }
        } else {
          double x = StringToDouble(tokens[0]);
          double y = StringToDouble(tokens[1]);
          switch (state) {
            case Polyline:
              out_polyline.push_back(p2t::Point(x, y));
              break;
            case Hole:
              assert(hole != nullptr);
              hole->push_back(p2t::Point(x, y));
              break;
            case Steiner:
              out_steiner.push_back(p2t::Point(x, y));
              break;
            default:
              assert(0);
          }
        }
      }
    } else {
      throw std::runtime_error("File not opened");
    }
  } catch (std::exception& e) {
    std::cerr << "Error parsing file: " << e.what() << std::endl;
    return false;
  }
  return true;
}

std::pair<p2t::Point, p2t::Point> BoundingBox(const std::vector<p2t::Point>& polyline)
{
  assert(polyline.size() > 0);
  using Scalar = decltype(p2t::Point::x);
  p2t::Point min(std::numeric_limits<Scalar>::max(), std::numeric_limits<Scalar>::max());
  p2t::Point max(std::numeric_limits<Scalar>::min(), std::numeric_limits<Scalar>::min());
  for (const p2t::Point& point : polyline) {
    min.x = std::min(min.x, point.x);
    min.y = std::min(min.y, point.y);
    max.x = std::max(max.x, point.x);
    max.y = std::max(max.y, point.y);
  }
  return std::make_pair(min, max);
}

void GenerateRandomPointDistribution(size_t num_points, double min, double max,
                                     std::vector<p2t::Point>& out_polyline,
                                     std::vector<std::vector<p2t::Point>>& out_holes, std::vector<p2t::Point>& out_steiner)
{
  out_polyline.push_back(p2t::Point(min, min));
  out_polyline.push_back(p2t::Point(min, max));
  out_polyline.push_back(p2t::Point(max, max));
  out_polyline.push_back(p2t::Point(max, min));

  max -= (1e-4);
  min += (1e-4);
  for (int i = 0; i < num_points; i++) {
    double x = Random(Fun, min, max);
    double y = Random(Fun, min, max);
    out_steiner.push_back(p2t::Point(x, y));
  }
}

void Init(int window_width, int window_height)
{
  if (glfwInit() != GL_TRUE)
    ShutDown(1);
  // width x height, 16 bit color, no depth, alpha or stencil buffers, windowed
  window = glfwCreateWindow(window_width, window_height, "Poly2Tri - C++", nullptr, nullptr);
  if (!window)
    ShutDown(1);

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
}

void ShutDown(int return_code)
{
  glfwTerminate();
  exit(return_code);
}

void MainLoop(double initial_zoom)
{
  // Zoom can be changed with the arrow keys
  double zoom = initial_zoom;

  bool running = true;
  while (running) {
    glfwPollEvents();
    glfwGetFramebufferSize(window, &window_width, &window_height);
    glViewport(0, 0, window_width, window_height);

    // Check if the ESCAPE key was pressed or the window was closed
    running = !glfwGetKey(window, GLFW_KEY_ESCAPE) && !glfwWindowShouldClose(window);

    // Press the UP and DOWN keys to zoom in/out. Press BACKSPACE to reset zoom.
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
      zoom *= 1.05;
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
      zoom *= 0.95;
    if (glfwGetKey(window, GLFW_KEY_BACKSPACE) == GLFW_PRESS)
      zoom = initial_zoom;

    // Draw the scene
    Draw(zoom);

    // Swap back and front buffers
    glfwSwapBuffers(window);
  }
}

void ResetZoom(double zoom, double cx, double cy, double width, double height, bool flag_flip_y)
{
  const double left =   -1.0 * width  / zoom;
  const double right =   1.0 * width  / zoom;
  const double bottom = -1.0 * height / zoom;
  const double top =     1.0 * height / zoom;
  const double flip_y = flag_flip_y ? -1.0 : 1.0;

  // Reset viewport
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  // Reset ortho view
  glOrtho(left, right, bottom, top, 1.0, -1.0);
  glTranslated(-cx, flip_y * -cy, 0.0);
  glScaled(1.0, flip_y, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glDisable(GL_DEPTH_TEST);
  glLoadIdentity();

  // Clear the screen
  glClear(GL_COLOR_BUFFER_BIT);
}

void Draw(const double zoom)
{
  // reset zoom
  p2t::Point center = p2t::Point(cx, cy);
  ResetZoom(zoom, center.x, center.y, (double)window_width, (double)window_height, flag_flip_y);

  for (int i = 0; i < triangles.size(); i++) {
    p2t::Triangle& t = *triangles[i];
    const p2t::Point& a = *t.GetPoint(0);
    const p2t::Point& b = *t.GetPoint(1);
    const p2t::Point& c = *t.GetPoint(2);

    // Red
    glColor3f(1, 0, 0);

    glBegin(GL_LINE_LOOP);
    glVertex2d(a.x, a.y);
    glVertex2d(b.x, b.y);
    glVertex2d(c.x, c.y);
    glEnd();
  }

  if (!convex_hull_triangulation)
  {
    // Green for constrained edges
    glColor3f(0, 1, 0);

    std::vector<std::vector<p2t::Point>*> polylines;
    polylines.push_back(&polyline);
    for (std::vector<p2t::Point>& hole : holes) {
        polylines.push_back(&hole);
    }
    for(int i = 0; i < polylines.size(); i++) {
      const std::vector<p2t::Point>& poly = *polylines[i];
      glBegin(GL_LINE_LOOP);
        for(int j = 0; j < poly.size(); j++) {
          glVertex2d(poly[j].x, poly[j].y);
        }
      glEnd();
    }
  }
}

void ConstrainedColor(bool constrained)
{
  if (constrained) {
    // Green
    glColor3f(0.f, 1.f, 0.f);
  } else {
    // Red
    glColor3f(1.f, 0.f, 0.f);
  }
}

double StringToDouble(const std::string& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return 0;
  return x;
}

double Fun(double x)
{
  return 2.5 + sin(10 * x) / x;
}

double Random(double (*fun)(double), double xmin = 0, double xmax = 1)
{
  static double (*Fun)(double) = nullptr, YMin, YMax;
  static bool First = true;

  // Initialises random generator for first call
  if (First)
  {
    First = false;
    srand((unsigned) time(nullptr));
  }

  // Evaluates maximum of function
  if (fun != Fun)
  {
    Fun = fun;
    YMin = 0, YMax = Fun(xmin);
    for (int iX = 1; iX < RAND_MAX; iX++)
    {
      double X = xmin + (xmax - xmin) * iX / RAND_MAX;
      double Y = Fun(X);
      YMax = Y > YMax ? Y : YMax;
    }
  }

  // Gets random values for X & Y
  double X = xmin + (xmax - xmin) * rand() / RAND_MAX;
  double Y = YMin + (YMax - YMin) * rand() / RAND_MAX;

  // Returns if valid and try again if not valid
  return Y < fun(X) ? X : Random(Fun, xmin, xmax);
}
