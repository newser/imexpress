#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <poly/cubic.h>

using namespace std;

TEST_CASE("poly_cubic")
{
    // 3 root
    auto x1 = iexp::poly::solve_cubic(1, -57.0, 1071.0, -6647.0);
    REQUIRE(__D_EQ_IN(get<0>(x1), 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<1>(x1), 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<2>(x1), 23.0, 1e-9));

    x1 = iexp::poly::solve_cubic(1, -109.0, 803.0, 50065.0);
    REQUIRE(__D_EQ_IN(get<0>(x1), -17.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<1>(x1), 31.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<2>(x1), 95.0, 1e-9));

    // 3 equal root
    x1 = iexp::poly::solve_cubic(1, -3, 3, -1);
    REQUIRE(__D_EQ_IN(get<0>(x1), 1, 1e-9));
    REQUIRE(__D_EQ_IN(get<0>(x1), 1, 1e-9));
    REQUIRE(__D_EQ_IN(get<0>(x1), 1, 1e-9));

    x1 = iexp::poly::solve_cubic(1, 0, 0, 1);
    REQUIRE(__D_EQ_IN(get<0>(x1), -1, 1e-9));
    REQUIRE(isnan(get<1>(x1)));
    REQUIRE(isnan(get<2>(x1)));
}

TEST_CASE("poly_complex_cubic")
{
    // 3 root
    auto x1 = iexp::poly::complex_solve_cubic(1, -57.0, 1071.0, -6647.0);
    REQUIRE(__D_EQ_IN(get<0>(x1).real(), 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<0>(x1).imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(get<1>(x1).real(), 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<1>(x1).imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(get<2>(x1).real(), 23.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<2>(x1).imag(), 0, 1e-9));

    x1 = iexp::poly::complex_solve_cubic(2, -2.0, 2.0, 78.0);
    REQUIRE(__D_EQ_IN(get<0>(x1).real(), -3.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<0>(x1).imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(get<1>(x1).real(), 2.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<1>(x1).imag(), -3.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<2>(x1).real(), 2.0, 1e-9));
    REQUIRE(__D_EQ_IN(get<2>(x1).imag(), 3.0, 1e-9));
}
