#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/wigner_symbol.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_wigner3j")
{
    iexp::Array<int, 6, 2> c;
    iexp::Array<double, 1, 2> rc, re;
    gsl_sf_result r;

    c.col(0) << 2, 4, 6, 0, 2, -2;
    c.col(1) << 4, 4, 8, 0, 0, 0;
    rc = sf::wigner3j(c);
    REQUIRE(__D_EQ9(rc(0), sqrt(8.0 / 105.0)));
    REQUIRE(__D_EQ9(rc(1), sqrt(2.0 / 35.0)));

    rc = sf::wigner3j(c, re);
    REQUIRE(__D_EQ9(rc(0), sqrt(8.0 / 105.0)));
    REQUIRE(__D_EQ9(rc(1), sqrt(2.0 / 35.0)));

    gsl_sf_coupling_3j_e(4, 4, 8, 0, 0, 0, &r);
    REQUIRE(__D_EQ9(re(1), r.err));
}

TEST_CASE("test_wigner6j")
{
    iexp::Array<int, 6, 1> c;
    iexp::Array<double, 1, 1> rc, re;
    gsl_sf_result r;

    c.col(0) << 2, 2, 4, 2, 2, 2;
    rc = sf::wigner6j(c);
    REQUIRE(__D_EQ9(rc(0), 1.0 / 6.0));

    rc = sf::wigner6j(c, re);
    REQUIRE(__D_EQ9(rc(0), 1.0 / 6.0));
    gsl_sf_coupling_6j_e(2, 2, 4, 2, 2, 2, &r);
    REQUIRE(__D_EQ9(re(0), r.err));
}

TEST_CASE("test_wigner9j")
{
    iexp::Array<int, 9, 1> c;
    iexp::Array<double, 1, 1> rc, re;
    gsl_sf_result r;

    c << 4, 2, 4, 3, 3, 2, 1, 1, 2;
    rc = sf::wigner9j(c);
    REQUIRE(__D_EQ9(rc(0), -sqrt(1.0 / 6.0) / 10.0));

    rc = sf::wigner9j(c, re);
    REQUIRE(__D_EQ9(rc(0), -sqrt(1.0 / 6.0) / 10.0));
    gsl_sf_coupling_9j_e(4, 2, 4, 3, 3, 2, 1, 1, 2, &r);
    REQUIRE(__D_EQ9(re(0), r.err));
}
