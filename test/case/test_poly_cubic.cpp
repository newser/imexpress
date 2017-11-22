#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <polynomial/cubic.h>

TEST_CASE("poly_cubic")
{
    iexp::ArrayXd c(3), x(3);

    c << -57.0, 1071.0, -6647.0;
    x = iexp::poly::solve_cubic(c);
    REQUIRE(__D_EQ_IN(x[0], 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(x[1], 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(x[2], 23.0, 1e-9));

    c << -109.0, 803.0, 50065.0;
    x = iexp::poly::solve_cubic(c);
    REQUIRE(__D_EQ_IN(x[0], -17.0, 1e-9));
    REQUIRE(__D_EQ_IN(x[1], 31.0, 1e-9));
    REQUIRE(__D_EQ_IN(x[2], 95.0, 1e-9));

    // compile
    x = iexp::poly::solve_cubic((c + c) * 2);
}

TEST_CASE("poly_complex_cubic")
{
    iexp::ArrayXd c(3);
    iexp::ArrayXcd x2(3);

    c << -57.0, 1071.0, -6647.0;
    x2 = iexp::poly::complex_solve_cubic(c);
    REQUIRE(__D_EQ_IN(x2[0].real(), 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[0].imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].real(), 17.0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[2].real(), 23.0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[2].imag(), 0, 1e-9));

    c << -1.0, 1.0, 39.0;
    x2 = iexp::poly::complex_solve_cubic(c);
    REQUIRE(__D_EQ_IN(x2[0].real(), -3.0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[0].imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].real(), 2.0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].imag(), -3.0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[2].real(), 2.0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[2].imag(), 3.0, 1e-9));

    // compile
    x2 = iexp::poly::complex_solve_cubic((c + c) * 2);
}
