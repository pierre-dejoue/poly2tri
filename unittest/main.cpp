#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Poly2triTest

#include "utility.h"

#include <poly2tri/poly2tri.h>

#include <boost/filesystem/path.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>

BOOST_AUTO_TEST_CASE(BasicTest)
{
  std::vector<p2t::Point*> polyline{
    new p2t::Point(0, 0),
    new p2t::Point(1, 0),
    new p2t::Point(1, 1),
  };
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 1);
  BOOST_REQUIRE(result[0] != nullptr);
  p2t::Triangle& t = *result[0];
  BOOST_CHECK_EQUAL(*t.GetPoint(0), *polyline[0]);
  BOOST_CHECK_EQUAL(*t.GetPoint(1), *polyline[1]);
  BOOST_CHECK_EQUAL(*t.GetPoint(2), *polyline[2]);
  BOOST_CHECK(t.IsInterior());
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(EdgeCases_EmptyInput)
{
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 0);
}

BOOST_AUTO_TEST_CASE(EdgeCases_SinglePoint)
{
  std::vector<p2t::Point*> polyline {
    new p2t::Point(1, 0)
  };
  p2t::CDT cdt;
  BOOST_CHECK_THROW(cdt.AddPolyline(polyline), std::invalid_argument);
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(EdgeCases_SingleEdge)
{
  std::vector<p2t::Point*> polyline {
    new p2t::Point(1, 0),
    new p2t::Point(0, 1)
  };
  p2t::CDT cdt;
  BOOST_CHECK_THROW(cdt.AddPolyline(polyline), std::invalid_argument);
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(EdgeCases_SingleSteinerPoint)
{
  p2t::Point steiner(1, 0);
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPoint(&steiner));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 0);
}

BOOST_AUTO_TEST_CASE(QuadTest)
{
  std::vector<p2t::Point*> polyline{ new p2t::Point(0, 0), new p2t::Point(0, 1),
                                     new p2t::Point(1, 1), new p2t::Point(1, 0) };
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPolyline(polyline));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(QuadSteinerTest)
{
  std::vector<p2t::Point*> points{ new p2t::Point(0, 0), new p2t::Point(0, 1),
                                   new p2t::Point(1, 1), new p2t::Point(1, 0) };
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPoints(points));
  BOOST_CHECK_NO_THROW(cdt.Triangulate(p2t::Policy::ConvexHull));
  const auto map = cdt.GetMap();
  BOOST_REQUIRE_EQUAL(map.size(), 6);
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  for (const auto p : points) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(QuadTestModernAPI)
{
  std::vector<p2t::Point> polyline{ p2t::Point(0, 0), p2t::Point(0, 1),
                                    p2t::Point(1, 1), p2t::Point(1, 0) };
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPolyline(polyline.data(), polyline.size()));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(QuadTestWithStride)
{
  struct PointStruct
  {
    p2t::Point p;
    float height;
    unsigned char type;
  };
  std::vector<PointStruct> polyline{
    PointStruct{ p2t::Point(0, 0), 0.1f, 245 },
    PointStruct{ p2t::Point(0, 1), 0.2f, 244 },
    PointStruct{ p2t::Point(1, 1), 0.3f, 244 },
    PointStruct{ p2t::Point(1, 0), 0.4f, 245 }
  };
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPolyline(reinterpret_cast<const p2t::Point*>(polyline.data()), polyline.size(), sizeof(PointStruct)));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(NarrowQuadTest)
{
  // Very narrow quad that used to demonstrate a failure case during
  // triangulation
  std::vector<p2t::Point*> polyline {
    new p2t::Point(0.0,     0.0),
    new p2t::Point(1.0e-05, 0.0),
    new p2t::Point(1.1e-04, 3.0e-07),
    new p2t::Point(1.0e-04, 3.0e-07)
  };
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(ConcaveBoundaryTest)
{
  // Concave-by-less-than-epsilon boundaries used to potentially fail
  // during triangulation
  const double eps = 1e-15; // This gave EdgeEvent - null triangle
  std::vector<p2t::Point*> polyline {
    new p2t::Point(0,0),
    new p2t::Point(0.5,eps),
    new p2t::Point(1,0),
    new p2t::Point(1-eps,0.836541),
    new p2t::Point(1,2),
    new p2t::Point(.46,1.46+eps),
    new p2t::Point(0,1),
    new p2t::Point(eps,0.5)
  };

  const double r2o4 = std::sqrt(2.)/4;
  std::vector<p2t::Point*> hole {
    new p2t::Point(0.5+r2o4,0.5),
    new p2t::Point(0.5,0.5+r2o4),
    new p2t::Point(0.5-r2o4,0.5),
    new p2t::Point(0.5-eps,0.5-r2o4)
  };

  std::vector<p2t::Point> interior_points
    {{0.21,0.79},{0.21,0.21},{0.79,0.21}};

  p2t::CDT cdt{ polyline };
  cdt.AddHole(hole);

  for (auto & p : interior_points)
    cdt.AddPoint(&p);

  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 18);
  BOOST_CHECK(IsConstrainedDelaunay(result));
  for (const auto p : polyline) {
    delete p;
  }
  for (const auto p : hole) {
    delete p;
  }
}


BOOST_AUTO_TEST_CASE(PolygonTest01)
{
  // Reported in issue #10
  std::vector<p2t::Point*> polyline = {
    new p2t::Point(-0.388419120000000006598384061363, 0.0368141516905975269002837535481),
    new p2t::Point(-0.388419120000000006598384061363, 0.0104235565411950892311665484158),
    new p2t::Point(-0.611580879999999993401615938637, 0.0104235565411950892311665484158),
    new p2t::Point(-0.611580879999999993401615938637, 0.1483950316905975341796875),
    new p2t::Point(-0.578899596898762469621146919962, 0.227294628589359948289683188705),
    new p2t::Point(-0.500000000000000000000000000000, 0.259975911690597527581303438637),
    new p2t::Point(+0.500000000000000000000000000000, 0.259975911690597527581303438637),
    new p2t::Point(+0.578899596898762469621146919962, 0.227294628589359948289683188705),
    new p2t::Point(+0.611580879999999993401615938637, 0.1483950316905975341796875),
    new p2t::Point(+0.611580879999999993401615938637, 0.0104235565411950614755909327869),
    new p2t::Point(+0.388419120000000006598384061363, 0.0104235565411950892311665484158),
    new p2t::Point(+0.388419120000000006598384061363, 0.0368141516905975130224959457337)
  };
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(result));
  // BOOST_CHECK(IsConstrainedDelaunay(result));            Fails
  BOOST_REQUIRE_EQUAL(result.size(), 10);
  for (const auto p : polyline) {
     delete p;
  }
}


BOOST_AUTO_TEST_CASE(PolygonTest02)
{
  // Reported in issue #10
  std::vector<p2t::Point*> polyline = {
    new p2t::Point(0.9636984967276516, 0.7676550649687783),
    new p2t::Point(0.9636984967276516, -0.7676550649687641),
    new p2t::Point(-0.3074475690811459, -0.7676550649687641),
    new p2t::Point(0.09401654924378076, -0.2590574983578904),
    new p2t::Point(0.10567230819363671, -0.09864698028880525),
    new p2t::Point(-0.03901177977841874, -0.028405214140875046),
    new p2t::Point(-0.428964921810446, -0.08483619470406722),
    new p2t::Point(-0.5128305980156834, -0.12847817634298053),
    new p2t::Point(-0.5512747518916774, -0.2148501697175078),
    new p2t::Point(-0.5917836778064418, -0.7037530067555622),
    new p2t::Point(-0.5520451065921502, -0.7676550649687641),
    new p2t::Point(-0.9636984967276516, -0.7676550649687641),
    new p2t::Point(-0.9636984967276516, 0.767655064968778)
  };
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_REQUIRE_EQUAL(result.size(), 11);
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(PolygonTest03)
{
  // Reported in issue #10
  std::vector<p2t::Point*> polyline = {
    new p2t::Point(0.9776422201600001, 0.9776422201599928),
    new p2t::Point(0.9776422201599999, -0.977642220160007),
    new p2t::Point(-0.12788518519240472, -0.9776422201599928),
    new p2t::Point(-0.3913394510746002, -0.33861494064331055),
    new p2t::Point(-0.47812835166211676, -0.9776422201599928),
    new p2t::Point(-0.9776422201600001, -0.9776422201599928),
    new p2t::Point(-0.9776422201600001, 0.977642220160007)
  };
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_REQUIRE_EQUAL(result.size(), 5);
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(PolygonTest04)
{
  std::vector<p2t::Point*> polyline {
    new p2t::Point(450, 2250),
    new p2t::Point(450, 1750),
    new p2t::Point(400, 1700),
    new p2t::Point(350, 1650),
    new p2t::Point(350, 500),
    new p2t::Point(1050, 1700)
  };

  std::vector<p2t::Point*> hole {
    new p2t::Point(980, 1636),
    new p2t::Point(950, 1600),
    new p2t::Point(650, 1230),
    new p2t::Point(625, 1247),
    new p2t::Point(600, 1250),
    new p2t::Point(591, 1350),
    new p2t::Point(550, 2050)
  };

  p2t::CDT cdt{ polyline };
  cdt.AddHole(hole);

  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_REQUIRE_EQUAL(result.size(), 13);
  for (const auto p : polyline) {
    delete p;
  }
  for (const auto p : hole) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(TestbedFilesTest)
{
  for (const auto& filename : { "custom.dat", "diamond.dat", "star.dat", "test.dat" }) {
    std::vector<p2t::Point*> polyline;
    // Load pointset from file
    // Parse and tokenize data file
    std::string line;
#ifndef P2T_BASE_DIR
    const auto basedir = boost::filesystem::path(__FILE__).remove_filename().parent_path();
#else
    const auto basedir = boost::filesystem::path(P2T_BASE_DIR);
#endif
    const auto datafile = basedir / boost::filesystem::path("testbed/data") / boost::filesystem::path(filename);
    std::ifstream myfile(datafile.string());
    BOOST_REQUIRE(myfile.is_open());
    while (!myfile.eof()) {
      getline(myfile, line);
      if (line.empty()) {
        break;
      }
      std::istringstream iss(line);
      std::vector<std::string> tokens;
      copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
           std::back_inserter<std::vector<std::string>>(tokens));
      double x = std::stod(tokens[0]);
      double y = std::stod(tokens[1]);
      polyline.push_back(new p2t::Point(x, y));
    }
    {
      const std::string case_message = std::string(filename) + " OuterPolygon " + std::to_string(polyline.size());
      p2t::CDT cdt{ polyline };
      BOOST_CHECK_NO_THROW(cdt.Triangulate(p2t::Policy::OuterPolygon));
      const auto result = cdt.GetTriangles();
      BOOST_REQUIRE(result.size() * 3 > polyline.size());
      BOOST_CHECK(TriangulationSanityChecks(result));
      BOOST_CHECK_MESSAGE(TriangulationSanityChecks(result), case_message);
      BOOST_CHECK_MESSAGE(IsConstrainedDelaunay(result), case_message);
    }
    {
      const std::string case_message = std::string(filename) + " ConvexHull " + std::to_string(polyline.size());
      p2t::CDT cdt{ polyline };
      BOOST_CHECK_NO_THROW(cdt.Triangulate(p2t::Policy::ConvexHull));
      const auto result = cdt.GetTriangles();
      BOOST_REQUIRE(result.size() * 3 > polyline.size());
      BOOST_CHECK(TriangulationSanityChecks(result));
      BOOST_CHECK_MESSAGE(TriangulationSanityChecks(result), case_message);
      BOOST_CHECK_MESSAGE(IsConstrainedDelaunay(result), case_message);
    }
    for (const auto p : polyline) {
      delete p;
    }
  }
}
