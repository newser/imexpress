#include <catch.hpp>
#include <iostream>
#include <math/inf_nan.h>

TEST_CASE("math_inf_nan")
{
    REQUIRE(IEXP_POSINF == GSL_POSINF);
    REQUIRE(IEXP_NEGINF == GSL_NEGINF);
    // REQUIRE(IEXP_NAN == GSL_NAN);

    SECTION("pos inf")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Zero();
        m(1, 1) = IEXP_POSINF;
        REQUIRE(iexp::isinf(m.array())(0, 0) == false);
        REQUIRE(iexp::isinf(m.array())(1, 1) == true);

        REQUIRE(iexp::isfinite(m.array())(0, 0) == true);
        REQUIRE(iexp::isfinite(m.array())(1, 1) == false);
    }

    SECTION("neg inf")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Zero();
        m(0, 1) = IEXP_NEGINF;
        REQUIRE(iexp::isinf(m.array())(0, 0) == false);
        REQUIRE(iexp::isinf(m.array())(0, 1) == true);

        REQUIRE(iexp::isfinite(m.array())(0, 0) == true);
        REQUIRE(iexp::isfinite(m.array())(0, 1) == false);
    }

    SECTION("NaN")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Zero();
        m(0, 1) = IEXP_NAN;
        REQUIRE(iexp::isnan(m.array())(0, 0) == false);
        REQUIRE(iexp::isnan(m.array())(0, 1) == true);
    }
}
