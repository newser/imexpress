#include <catch.hpp>
#include <fmin/fmin.h>
#include <fmin/mfmin.h>
#include <iostream>
#include <math/constant.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_fmin")
{
    iexp::fmin m([](double x, void *) { return pow(x, 4) - 1; }, -3, 17, -1);
    double r = m.find(0.001, 0.001);
    REQUIRE(__D_EQ3(r, 0));

    iexp::fmin m2([](double x, void *) { return cos(x); }, 0, 6, 3);
    r = m2.find(0.001, 0.001);
    REQUIRE(__D_EQ3(r, M_PI));
}
