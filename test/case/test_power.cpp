#include <catch.hpp>
#include <iostream>
#include <special/power.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_power")
{
    Matrix2d p, r, e;
    p << 1, 2, 3, 4;
    r = sf::pow(2, p);
    REQUIRE(r(0, 0) == 1 * 1);
    REQUIRE(r(0, 1) == 2 * 2);
    REQUIRE(r(1, 0) == 3 * 3);
    REQUIRE(r(1, 1) == 4 * 4);

    r = sf::pow(2, p, e);
    REQUIRE(r(0, 0) == 1 * 1);
    REQUIRE(r(0, 1) == 2 * 2);
    REQUIRE(r(1, 0) == 3 * 3);
    REQUIRE(r(1, 1) == 4 * 4);
}
