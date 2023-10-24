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
  std::vector<p2t::Point> points {
    p2t::Point(0, 0),
    p2t::Point(1, 0),
    p2t::Point(1, 1),
  };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 1);
  BOOST_REQUIRE(result[0] != nullptr);
  p2t::Triangle& t = *result[0];
  BOOST_CHECK_EQUAL(*t.GetPoint(0), *polyline[0]);
  BOOST_CHECK_EQUAL(*t.GetPoint(1), *polyline[1]);
  BOOST_CHECK_EQUAL(*t.GetPoint(2), *polyline[2]);
  BOOST_CHECK(t.IsInterior());
}

BOOST_AUTO_TEST_CASE(EdgeCases_EmptyInput)
{
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 0);
}

BOOST_AUTO_TEST_CASE(EdgeCases_SinglePoint)
{
  std::vector<p2t::Point> points {
    p2t::Point(1, 0)
  };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  BOOST_CHECK_THROW(cdt.AddPolyline(polyline), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(EdgeCases_SingleEdge)
{
  std::vector<p2t::Point> points {
    p2t::Point(1, 0),
    p2t::Point(0, 1)
  };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  BOOST_CHECK_THROW(cdt.AddPolyline(polyline), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(EdgeCases_SingleSteinerPoint)
{
  p2t::Point steiner(1, 0);
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPoint(&steiner));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 0);
}

BOOST_AUTO_TEST_CASE(QuadTest)
{
  std::vector<p2t::Point> points { p2t::Point(0, 0), p2t::Point(0, 1),
                                   p2t::Point(1, 1), p2t::Point(1, 0) };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPolyline(polyline));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(QuadSteinerTest)
{
  std::vector<p2t::Point> points { p2t::Point(0, 0), p2t::Point(0, 1),
                                   p2t::Point(1, 1), p2t::Point(1, 0) };
  const auto steiner_points = MakePointerVector(points);
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPoints(steiner_points));
  BOOST_CHECK_NO_THROW(cdt.Triangulate(p2t::Policy::ConvexHull));
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(QuadTestModernAPI)
{
  std::vector<p2t::Point> polyline{ p2t::Point(0, 0), p2t::Point(0, 1),
                                    p2t::Point(1, 1), p2t::Point(1, 0) };
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPolyline(polyline.data(), polyline.size()));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
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
  std::vector<PointStruct> polyline {
    PointStruct{ p2t::Point(0, 0), 0.1f, 245 },
    PointStruct{ p2t::Point(0, 1), 0.2f, 244 },
    PointStruct{ p2t::Point(1, 1), 0.3f, 244 },
    PointStruct{ p2t::Point(1, 0), 0.4f, 245 }
  };
  p2t::CDT cdt;
  BOOST_CHECK_NO_THROW(cdt.AddPolyline(reinterpret_cast<const p2t::Point*>(polyline.data()), polyline.size(), sizeof(PointStruct)));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(MultipleTriangulations)
{
  std::vector<p2t::Point> polyline{ p2t::Point(0, 0), p2t::Point(0, 3),
                                    p2t::Point(3, 3), p2t::Point(3, 0) };
  std::vector<p2t::Point> points  { p2t::Point(1, 1), p2t::Point(1, 2),
                                    p2t::Point(2, 2), p2t::Point(2, 1) };

  // Allocate one CDT object
  p2t::CDT cdt;

  // First triangulation
  BOOST_CHECK_NO_THROW(cdt.AddPolyline(polyline.data(), polyline.size()));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  auto result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);

  // Add Steiner and perform another triangulation
  BOOST_CHECK_NO_THROW(cdt.AddPoints(points.data(), points.size()));
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 10);
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(NarrowQuadTest)
{
  // Very narrow quad that used to demonstrate a failure case during
  // triangulation
  std::vector<p2t::Point> points {
    p2t::Point(0.0,     0.0),
    p2t::Point(1.0e-05, 0.0),
    p2t::Point(1.1e-04, 3.0e-07),
    p2t::Point(1.0e-04, 3.0e-07)
  };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 2);
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(Pyramide)
{
  // Test case with PointEvents where the new point's projection on the advancing front coincides exactly with a node
  std::vector<p2t::Point> outer_points {
    p2t::Point(0.0, 0.0),
    p2t::Point(2.0, 0.0),
    p2t::Point(1.0, 3.0)
  };
  std::vector<p2t::Point> points {
    p2t::Point(1.0, 1.0),
    p2t::Point(1.0, 2.0)
  };
  p2t::CDT cdt;
  cdt.AddPolyline(outer_points.data(), outer_points.size());
  cdt.AddPoints(points.data(), points.size());
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(result.size(), 5);
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
}

BOOST_AUTO_TEST_CASE(ConcaveBoundaryTest)
{
  // Concave-by-less-than-epsilon boundaries used to potentially fail
  // during triangulation
  const double eps = 1e-15; // This gave EdgeEvent - null triangle
  std::vector<p2t::Point> points {
    p2t::Point(0,0),
    p2t::Point(0.5,eps),
    p2t::Point(1,0),
    p2t::Point(1-eps,0.836541),
    p2t::Point(1,2),
    p2t::Point(.46,1.46+eps),
    p2t::Point(0,1),
    p2t::Point(eps,0.5)
  };

  const double r2o4 = std::sqrt(2.)/4;
  std::vector<p2t::Point> hole_points {
    p2t::Point(0.5+r2o4,0.5),
    p2t::Point(0.5,0.5+r2o4),
    p2t::Point(0.5-r2o4,0.5),
    p2t::Point(0.5-eps,0.5-r2o4)
  };

  std::vector<p2t::Point> interior_points {
    p2t::Point(0.21,0.79),
    p2t::Point(0.21,0.21),
    p2t::Point(0.79,0.21)
  };

  const auto polyline = MakePointerVector(points);
  const auto hole = MakePointerVector(hole_points);

  p2t::CDT cdt;
  cdt.AddPolyline(polyline);
  cdt.AddHole(hole);

  for (auto & p : interior_points)
    cdt.AddPoint(&p);

  BOOST_CHECK_THROW(cdt.Triangulate(), std::runtime_error);
}

namespace assets {
  std::vector<p2t::Point>& polygon_test_01()
  {
    // Reported in issue #10
    static std::vector<p2t::Point> points {
      p2t::Point(-0.388419120000000006598384061363, 0.0368141516905975269002837535481),
      p2t::Point(-0.388419120000000006598384061363, 0.0104235565411950892311665484158),
      p2t::Point(-0.611580879999999993401615938637, 0.0104235565411950892311665484158),
      p2t::Point(-0.611580879999999993401615938637, 0.1483950316905975341796875),
      p2t::Point(-0.578899596898762469621146919962, 0.227294628589359948289683188705),
      p2t::Point(-0.500000000000000000000000000000, 0.259975911690597527581303438637),
      p2t::Point(+0.500000000000000000000000000000, 0.259975911690597527581303438637),
      p2t::Point(+0.578899596898762469621146919962, 0.227294628589359948289683188705),
      p2t::Point(+0.611580879999999993401615938637, 0.1483950316905975341796875),
      p2t::Point(+0.611580879999999993401615938637, 0.0104235565411950892311665484158),
      p2t::Point(+0.388419120000000006598384061363, 0.0104235565411950892311665484158),
      p2t::Point(+0.388419120000000006598384061363, 0.0368141516905975130224959457337)
    };
    return points;
  }
} // namespace assets

BOOST_AUTO_TEST_CASE(PolygonTest01)
{
  auto& points = assets::polygon_test_01();
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  cdt.AddPolyline(polyline);
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  // BOOST_CHECK(IsConstrainedDelaunay(result));            Fails
  BOOST_REQUIRE_EQUAL(result.size(), 10);
}

BOOST_AUTO_TEST_CASE(PolygonTest01_ConvexHull)
{
  auto& points = assets::polygon_test_01();
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  cdt.AddPolyline(polyline);
  BOOST_CHECK_NO_THROW(cdt.Triangulate(p2t::Policy::ConvexHull));
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  // BOOST_CHECK(IsConstrainedDelaunay(result));            Fails
  BOOST_REQUIRE_EQUAL(result.size(), 12);
}

BOOST_AUTO_TEST_CASE(PolygonTest02)
{
  // Reported in issue #10
  std::vector<p2t::Point> points {
    p2t::Point(0.9636984967276516, 0.7676550649687783),
    p2t::Point(0.9636984967276516, -0.7676550649687641),
    p2t::Point(-0.3074475690811459, -0.7676550649687641),
    p2t::Point(0.09401654924378076, -0.2590574983578904),
    p2t::Point(0.10567230819363671, -0.09864698028880525),
    p2t::Point(-0.03901177977841874, -0.028405214140875046),
    p2t::Point(-0.428964921810446, -0.08483619470406722),
    p2t::Point(-0.5128305980156834, -0.12847817634298053),
    p2t::Point(-0.5512747518916774, -0.2148501697175078),
    p2t::Point(-0.5917836778064418, -0.7037530067555622),
    p2t::Point(-0.5520451065921502, -0.7676550649687641),
    p2t::Point(-0.9636984967276516, -0.7676550649687641),
    p2t::Point(-0.9636984967276516, 0.767655064968778)
  };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_REQUIRE_EQUAL(result.size(), 11);
}

BOOST_AUTO_TEST_CASE(PolygonTest03)
{
  // Reported in issue #10
  std::vector<p2t::Point> points {
    p2t::Point(0.9776422201600001, 0.9776422201599928),
    p2t::Point(0.9776422201599999, -0.977642220160007),
    p2t::Point(-0.12788518519240472, -0.9776422201599928),
    p2t::Point(-0.3913394510746002, -0.33861494064331055),
    p2t::Point(-0.47812835166211676, -0.9776422201599928),
    p2t::Point(-0.9776422201600001, -0.9776422201599928),
    p2t::Point(-0.9776422201600001, 0.977642220160007)
  };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt{ polyline };
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_REQUIRE_EQUAL(result.size(), 5);
}

BOOST_AUTO_TEST_CASE(PolygonTest04)
{
  std::vector<p2t::Point> points {
    p2t::Point(450, 2250),
    p2t::Point(450, 1750),
    p2t::Point(400, 1700),
    p2t::Point(350, 1650),
    p2t::Point(350, 500),
    p2t::Point(1050, 1700)
  };

  std::vector<p2t::Point> hole_points {
    p2t::Point(980, 1636),
    p2t::Point(950, 1600),
    p2t::Point(650, 1230),
    p2t::Point(625, 1247),
    p2t::Point(600, 1250),
    p2t::Point(591, 1350),
    p2t::Point(550, 2050)
  };

  const auto polyline = MakePointerVector(points);
  const auto hole = MakePointerVector(hole_points);
  p2t::CDT cdt{ polyline };
  cdt.AddHole(hole);

  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_REQUIRE_EQUAL(result.size(), 13);
}

BOOST_AUTO_TEST_CASE(PolygonTest05)
{
  std::vector<p2t::Point> points {
    p2t::Point(-0.37557949531446866,  -38.050538121782154),
    p2t::Point( 50.922527881862266,   -55.434201094814995),
    p2t::Point( 103.925537109375,     -48.979888916015625),
    p2t::Point( 145.67479299050069,   -20.598233547940595),
    p2t::Point( 166.02706909179688,    24.920921325683594),
    p2t::Point( 164.7850341796875,     98.821723937988281),
    p2t::Point( 27.873260498046875,    175.64570617675781),
  };
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  cdt.AddPolyline(polyline);
  BOOST_CHECK_NO_THROW(cdt.Triangulate());
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_REQUIRE_EQUAL(result.size(), 5);
}

namespace assets {
  std::vector<p2t::Point>& wavy_line()
  {
    static std::vector<p2t::Point> points {
      p2t::Point(401.15995152590949, 2.5053475482118301),
      p2t::Point(410.16662617780162, 17.270662776009431),
      p2t::Point(419.17332095751101, 32.036002863748990),
      p2t::Point(428.17999267578125, 46.801330566406250),
      p2t::Point(436.84804989727274, 61.263621915466842),
      p2t::Point(445.02004811036909, 76.010577151077314),
      p2t::Point(452.25342660665638, 91.236749191325998),
      p2t::Point(458.01225656235124, 107.07472583751746),
      p2t::Point(461.64973074721024, 123.52260505815084),
      p2t::Point(462.50531005859375, 140.33868408203125),
      p2t::Point(460.11581070403435, 157.70633963277308),
      p2t::Point(454.99565790092311, 174.48500610671999),
      p2t::Point(448.07518559972414, 190.61363696800072),
      p2t::Point(440.18489398711336, 206.29529570498016),
      p2t::Point(432.03607177734375, 221.84536743164062),
      p2t::Point(424.42463680717765, 237.65366111768546),
      p2t::Point(418.07214153119912, 254.00409074140316),
      p2t::Point(413.69684453908076, 270.98256618320261),
      p2t::Point(412.24510344656170, 288.44153751137412),
      p2t::Point(414.70538330078125, 305.77069091796875),
      p2t::Point(420.46156446287364, 320.59940607346948),
      p2t::Point(428.22786234750356, 334.50267741057746),
      p2t::Point(436.63321511421054, 348.03577234163731),
      p2t::Point(444.46611253001652, 361.90211653439542),
      p2t::Point(450.24195952386771, 376.72003964358720),
      p2t::Point(452.06939697265625, 392.48269653320312),
      p2t::Point(448.27310352469021, 409.42717564103714),
      p2t::Point(439.65798952164755, 424.55024261857318),
      p2t::Point(428.01873779296875, 437.52001953125000),
      p2t::Point(414.62341937267684, 448.70653659891650),
      p2t::Point(400.57632614128829, 459.07455910582212),
      p2t::Point(386.73602294921875, 469.71334838867188),
      p2t::Point(374.56233041482585, 480.81085246114424),
      p2t::Point(363.48344964792227, 493.00179063955875),
      p2t::Point(353.58105447115668, 506.16639598127131),
      p2t::Point(344.93440453844391, 520.18777249177572),
      p2t::Point(337.62332807852556, 534.94953152461972),
      p2t::Point(331.72930908203125, 550.33178710937500),
      p2t::Point(327.64024877079913, 564.27598675533295),
      p2t::Point(324.02377906389802, 578.35220682341105),
      p2t::Point(320.20061155933956, 592.37288271372699),
      p2t::Point(315.49603271484375, 606.11883544921875),
      p2t::Point(309.06802817859591, 620.10105318281865),
      p2t::Point(301.52470812013962, 633.51915309171090),
      p2t::Point(293.41170990285127, 646.60315997383395),
      p2t::Point(285.23134054601928, 659.64542999720902),
      p2t::Point(277.47735595703125, 672.94378662109375),
      p2t::Point(270.74953478199785, 686.77273580487781),
      p2t::Point(265.55216795161709, 701.24273774196627),
      p2t::Point(262.40670276470973, 716.28597580680525),
      p2t::Point(262.01081370798067, 731.63863848110520),
      p2t::Point(265.05865478515625, 746.67871093750000)
    };
    return points;
  }
} // namespace assets

BOOST_AUTO_TEST_CASE(AddOpenPolyline_ConvexHull)
{
  auto& points = assets::wavy_line();
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  cdt.AddOpenPolyline(polyline.data(), polyline.size());
  BOOST_CHECK_NO_THROW(cdt.Triangulate(p2t::Policy::ConvexHull));
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
  BOOST_CHECK(IsConstrainedDelaunay(result));
  BOOST_CHECK_EQUAL(result.size(), 88);
}

BOOST_AUTO_TEST_CASE(AddOpenPolyline_OuterPolygon)
{
  auto& points = assets::wavy_line();
  const auto polyline = MakePointerVector(points);
  p2t::CDT cdt;
  cdt.AddOpenPolyline(polyline.data(), polyline.size());
  BOOST_CHECK_NO_THROW(cdt.Triangulate(p2t::Policy::OuterPolygon));
  const auto& result = cdt.GetTriangles();
  BOOST_CHECK(TriangulationSanityChecks(cdt, result));
}

struct FileTest
{
  const char* filename;
  const p2t::Policy triangulation_policy;
  const bool expected_is_cdt;
  const std::size_t expected_points;
  const std::size_t expected_triangles;
};

BOOST_AUTO_TEST_CASE(TestbedFilesTest)
{
  const std::vector<FileTest> file_tests {
    { "diamond.dat", p2t::Policy::OuterPolygon, true, 10, 8 },
    { "diamond.dat", p2t::Policy::ConvexHull, true, 10, 8 },
    { "star.dat", p2t::Policy::OuterPolygon, true, 10, 8 },
    { "star.dat", p2t::Policy::ConvexHull, true, 10, 13 },
    { "test.dat", p2t::Policy::OuterPolygon, true, 6, 4 },
    { "test.dat", p2t::Policy::ConvexHull, true, 6, 6 },
    { "dude_sampled_1.dat", p2t::Policy::OuterPolygon, true, 109, 107 },
    { "dude_sampled_1.dat", p2t::Policy::ConvexHull, true, 109, 199 },
    { "dude_sampled_2.dat", p2t::Policy::OuterPolygon, true, 512, 510 },
    { "dude_sampled_2.dat", p2t::Policy::ConvexHull, true, 512, 998 },
  };
  for (const auto& test_case : file_tests) {
    std::vector<p2t::Point> points;
    // Load pointset from file
    // Parse and tokenize data file
    std::string line;
#ifndef P2T_BASE_DIR
    const auto basedir = boost::filesystem::path(__FILE__).remove_filename().parent_path();
#else
    const auto basedir = boost::filesystem::path(P2T_BASE_DIR);
#endif
    const auto datafile = basedir / boost::filesystem::path("testbed/data") / boost::filesystem::path(test_case.filename);
    std::ifstream myfile(datafile.string());
    BOOST_REQUIRE(myfile.is_open());
    while (!myfile.eof()) {
      std::getline(myfile, line);
      if (line.empty()) {
        break;
      }
      std::istringstream iss(line);
      std::vector<std::string> tokens;
      std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
           std::back_inserter<std::vector<std::string>>(tokens));
      double x = std::stod(tokens[0]);
      double y = std::stod(tokens[1]);
      points.push_back(p2t::Point(x, y));
    }
    const auto polyline = MakePointerVector(points);
    const std::string case_message = [&]() {
      std::stringstream sout;
      sout << test_case.filename << " " << test_case.triangulation_policy << " " << polyline.size();
      return sout.str();
    }();
    p2t::CDT cdt;
    cdt.AddPolyline(polyline);
    BOOST_CHECK_NO_THROW(cdt.Triangulate(test_case.triangulation_policy));
    const auto& result = cdt.GetTriangles();
    BOOST_REQUIRE(result.size() * 3 > polyline.size());
    BOOST_CHECK_MESSAGE(TriangulationSanityChecks(cdt, result), case_message);
    BOOST_CHECK_MESSAGE(IsConstrainedDelaunay(result) == test_case.expected_is_cdt, case_message);
    BOOST_CHECK_MESSAGE(result.size() == test_case.expected_triangles, case_message);
    BOOST_CHECK_MESSAGE(CountPoints(result) == test_case.expected_points, case_message);
  }
}
