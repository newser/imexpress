#include <catch.hpp>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math/constant.h>

TEST_CASE("error")
{
    bool ok = false;

    try {
        GSL_ERROR("test", GSL_ERANGE);
    } catch (std::exception &e) {
        if (strncmp(e.what(),
                    "output range error:test",
                    sizeof("output range error:test") - 1) == 0) {
            ok = true;
        }
    }
    REQUIRE(ok);
}

TEST_CASE("math_constant")
{
    REQUIRE(IEXP_E == M_E);
    REQUIRE(IEXP_LOG2E == M_LOG2E);
    REQUIRE(IEXP_LOG10E == M_LOG10E);
    REQUIRE(IEXP_SQRT2 == M_SQRT2);
    REQUIRE(IEXP_SQRT1_2 == M_SQRT1_2);
    REQUIRE(IEXP_SQRT3 == M_SQRT3);
    REQUIRE(IEXP_PI == M_PI);
    REQUIRE(IEXP_PI_2 == M_PI_2);
    REQUIRE(IEXP_PI_4 == M_PI_4);
    REQUIRE(IEXP_SQRTPI == M_SQRTPI);
    REQUIRE(IEXP_2_SQRTPI == M_2_SQRTPI);
    REQUIRE(IEXP_1_PI == M_1_PI);
    REQUIRE(IEXP_2_PI == M_2_PI);
    REQUIRE(IEXP_LN10 == M_LN10);
    REQUIRE(IEXP_LN2 == M_LN2);
    REQUIRE(IEXP_LNPI == M_LNPI);
    REQUIRE(IEXP_EULER == M_EULER);
}
