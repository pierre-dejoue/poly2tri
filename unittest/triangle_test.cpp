#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif
#include <boost/test/unit_test.hpp>
#include <poly2tri/common/shapes.h>

BOOST_AUTO_TEST_CASE(TriangleTest)
{
  p2t::Point a(0.0, 0.0);
  p2t::Point b(1.0, 0.0);
  p2t::Point c(0.5, 0.5);
  p2t::Point b_cpy = b;
  p2t::Triangle triangle(&a, &b, &c);
  BOOST_CHECK(triangle.Contains(&a));
  BOOST_CHECK(triangle.Contains(&b));
  BOOST_CHECK(triangle.Contains(&c));
  BOOST_CHECK(!triangle.Contains(&b_cpy));
  BOOST_CHECK(triangle.CircumcircleContains(p2t::Point(0.5, 0.1)));
  BOOST_CHECK(!triangle.CircumcircleContains(p2t::Point(1.0, 0.4)));
  BOOST_CHECK_EQUAL(triangle.GetOrientation(), p2t::Orientation::CCW);
}
