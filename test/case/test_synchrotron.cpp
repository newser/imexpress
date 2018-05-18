#include <catch.hpp>
#include <iostream>
#include <special/synchrotron.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_synchrotron")
{
    gsl_sf_result re;
    Matrix2d p, r, e;
    p << 0.01, 1, 10, 100;

    // synchrotron1
    r = sf::synchrotron1(p);
    REQUIRE(r(0, 0) == gsl_sf_synchrotron_1(p(0, 0)));
    REQUIRE(r(1, 1) == gsl_sf_synchrotron_1(p(1, 1)));

    r = sf::synchrotron1(p, e);
    gsl_sf_synchrotron_1_e(p(0, 0), &re);
    REQUIRE(r(0, 0) == re.val);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_synchrotron_1_e(p(1, 1), &re);
    REQUIRE(r(1, 1) == re.val);
    REQUIRE(re.err == e(1, 1));

    // synchrotron2
    r = sf::synchrotron2(p);
    REQUIRE(r(0, 0) == gsl_sf_synchrotron_2(p(0, 0)));
    REQUIRE(r(1, 1) == gsl_sf_synchrotron_2(p(1, 1)));

    r = sf::synchrotron2(p, e);
    gsl_sf_synchrotron_2_e(p(0, 0), &re);
    REQUIRE(r(0, 0) == re.val);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_synchrotron_2_e(p(1, 1), &re);
    REQUIRE(r(1, 1) == re.val);
    REQUIRE(re.err == e(1, 1));
}
