#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <poly/quadratic.h>

TEST_CASE("poly_quadratic")
{
    iexp::ArrayXd c(3), x(2);
    iexp::ArrayXcd x2(2);

    c << 4.0, -20.0, 26.0;
    x = iexp::poly::solve_quad(c);
    REQUIRE((iexp::isnan(x)[0] && iexp::isnan(x)[1]));

    x2 = iexp::poly::complex_solve_quad(c);
    REQUIRE(__D_EQ_IN(x2[0].real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(x2[0].imag(), -0.5, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].imag(), 0.5, 1e-9));

    c << 4.0, -20.0, 25.0;
    x = iexp::poly::solve_quad(c);
    REQUIRE(__D_EQ_IN(x[0], 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(x[1], 2.5, 1e-9));

    x2 = iexp::poly::complex_solve_quad(c);
    REQUIRE(__D_EQ_IN(x2[0].real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(x2[0].imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].imag(), 0, 1e-9));

    c << 4.0, -20.0, 21.0;
    x = iexp::poly::solve_quad(c);
    REQUIRE(__D_EQ_IN(x[0], 1.5, 1e-9));
    REQUIRE(__D_EQ_IN(x[1], 3.5, 1e-9));

    x2 = iexp::poly::complex_solve_quad(c);
    REQUIRE(__D_EQ_IN(x2[0].real(), 1.5, 1e-9));
    REQUIRE(__D_EQ_IN(x2[0].imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].real(), 3.5, 1e-9));
    REQUIRE(__D_EQ_IN(x2[1].imag(), 0, 1e-9));

    // compile
    x = iexp::poly::solve_quad((c + c) * 2);
    x2 = iexp::poly::complex_solve_quad((c + c) * 2);
}
