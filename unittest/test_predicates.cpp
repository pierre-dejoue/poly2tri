#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE PredicatesTest

#include <boost/test/unit_test.hpp>

#include <common/predicates.h>

template <typename FP = double>
struct Point {
  using scalar = FP;
  scalar x, y;

  Point(scalar x, scalar y) : x(x), y(y) {}
};

BOOST_AUTO_TEST_CASE(Orient2d_Basic)
{
  using namespace p2t;
  BOOST_CHECK(Orient2d(Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, -1.0)) == Orientation::CW);
  BOOST_CHECK(Orient2d(Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0,  0.0)) == Orientation::COLLINEAR);
  BOOST_CHECK(Orient2d(Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0,  1.0)) == Orientation::CCW);
}

BOOST_AUTO_TEST_CASE(InCircle_Basic)
{
  Point a(0.0, 0.0);
  Point b(1.0, 0.0);
  Point c(0.0, 2.0);
  Point d_inside(1.0, 1.0);
  Point d_border(1.0, 2.0);
  Point d_outside(1.0, 3.0);

  BOOST_CHECK(p2t::InCircle(a, b, c, d_inside) == true);
  BOOST_CHECK(p2t::InCircle(a, b, c, d_border) == false);
  BOOST_CHECK(p2t::InCircle(a, b, c, d_outside) == false);
}

BOOST_AUTO_TEST_CASE(InCircle_UseCase1)
{
  Point a(251.33397041666666, 526.43265294642856);
  Point b(251.33290750000000, 526.52353303571419);
  Point c(251.32971874999998, 526.79617330357132);
  Point d(251.33078166666667, 526.70529321428558);
  BOOST_CHECK(p2t::InCircle(a, b, c, d) == false);
  BOOST_CHECK(p2t::InCircle(c, a, d, b) == false);
}

BOOST_AUTO_TEST_CASE(InCircle_UseCase2)
{
  Point a(0.50000000000000000, 0.85355339059327373);
  Point b(0.14644660940672621, 0.50000000000000000);
  Point c(0.85355339059327373, 0.50000000000000000);
  Point d(0.49999999999999900, 0.14644660940672621);
  BOOST_CHECK(p2t::InCircle(a, b, c, d) == false);
  BOOST_CHECK(p2t::InCircle(c, a, d, b) == false);
}
