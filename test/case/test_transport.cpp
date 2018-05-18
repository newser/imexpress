#include <catch.hpp>
#include <iostream>
#include <special/transport.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_transport")
{
    gsl_sf_result re;
    Matrix2d p, r, e;
    p << 1.0e-10, 0.1, 3, 100;

    // 2
    r = sf::transport2(p);
    REQUIRE(r(0, 0) == gsl_sf_transport_2(p(0, 0)));
    REQUIRE(r(1, 1) == gsl_sf_transport_2(p(1, 1)));

    r = sf::transport2(p, e);
    gsl_sf_transport_2_e(p(0, 0), &re);
    REQUIRE(r(0, 0) == re.val);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_transport_2_e(p(1, 1), &re);
    REQUIRE(r(1, 1) == re.val);
    REQUIRE(re.err == e(1, 1));

    // 3
    r = sf::transport3(p);
    REQUIRE(r(0, 0) == gsl_sf_transport_3(p(0, 0)));
    REQUIRE(r(1, 1) == gsl_sf_transport_3(p(1, 1)));

    r = sf::transport3(p, e);
    gsl_sf_transport_3_e(p(0, 0), &re);
    REQUIRE(r(0, 0) == re.val);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_transport_3_e(p(1, 1), &re);
    REQUIRE(r(1, 1) == re.val);
    REQUIRE(re.err == e(1, 1));

    // 2
    r = sf::transport4(p);
    REQUIRE(r(0, 0) == gsl_sf_transport_4(p(0, 0)));
    REQUIRE(r(1, 1) == gsl_sf_transport_4(p(1, 1)));

    r = sf::transport4(p, e);
    gsl_sf_transport_4_e(p(0, 0), &re);
    REQUIRE(r(0, 0) == re.val);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_transport_4_e(p(1, 1), &re);
    REQUIRE(r(1, 1) == re.val);
    REQUIRE(re.err == e(1, 1));

    // 2
    r = sf::transport5(p);
    REQUIRE(r(0, 0) == gsl_sf_transport_5(p(0, 0)));
    REQUIRE(r(1, 1) == gsl_sf_transport_5(p(1, 1)));

    r = sf::transport5(p, e);
    gsl_sf_transport_5_e(p(0, 0), &re);
    REQUIRE(r(0, 0) == re.val);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_transport_5_e(p(1, 1), &re);
    REQUIRE(r(1, 1) == re.val);
    REQUIRE(re.err == e(1, 1));
}
