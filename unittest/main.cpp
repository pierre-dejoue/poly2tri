#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Poly2triTest

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
  BOOST_CHECK_EQUAL(*result[0]->GetPoint(0), *polyline[0]);
  BOOST_CHECK_EQUAL(*result[0]->GetPoint(1), *polyline[1]);
  BOOST_CHECK_EQUAL(*result[0]->GetPoint(2), *polyline[2]);
  BOOST_CHECK(p2t::IsDelaunay(result));
  for (const auto p : polyline) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(QuadTest)
{
  std::vector<p2t::Point*> polyline{ new p2t::Point(0, 0), new p2t::Point(0, 1),
                                     new p2t::Point(1, 1), new p2t::Point(1, 0) };
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(p2t::IsDelaunay(result));
  for (const auto p : polyline) {
    delete p;
  }
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
  BOOST_CHECK(p2t::IsDelaunay(result));
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
  BOOST_CHECK(p2t::IsDelaunay(result));
  for (const auto p : polyline) {
    delete p;
  }
  for (const auto p : hole) {
    delete p;
  }
}

BOOST_AUTO_TEST_CASE(SimpleTestCaseInfiniteLoop)
{
  std::vector<p2t::Point*> polyline {
    new p2t::Point(450, 2250),
    new p2t::Point(450, 1750),
    new p2t::Point(400, 1700),
    new p2t::Point(350, 1650),
    new p2t::Point(350, 500),
    new p2t::Point(1050, 1700)
  };

  std::vector<p2t::Point*> hole = {
    new p2t::Point(980, 1636),
    new p2t::Point(950, 1600),
    new p2t::Point(650, 1230),
    new p2t::Point(625, 1247),
    new p2t::Point(600, 1250),
    new p2t::Point(591, 1350),
    new p2t::Point(550, 2050)
  };

  p2t::CDT cdt(polyline);
  cdt.AddHole(hole);
  cdt.Triangulate();
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
    p2t::CDT cdt{ polyline };
    BOOST_CHECK_NO_THROW(cdt.Triangulate());
    const auto result = cdt.GetTriangles();
    BOOST_REQUIRE(result.size() * 3 > polyline.size());
    BOOST_CHECK_MESSAGE(p2t::IsDelaunay(result), std::string(filename) + " " + std::to_string(polyline.size()));
    for (const auto p : polyline) {
      delete p;
    }
  }
}
