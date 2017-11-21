#include <../test/test_util.h>
#include <Eigen/Dense>
#include <catch.hpp>
#include <iostream>
#include <polynomial/general.h>

TEST_CASE("poly_general")
{
    iexp::ArrayXd c(6);
    iexp::ArrayXcd x(5);

    c << -120, 274, -225, 85, -15, 1;
    x = iexp::poly::complex_solve(c);
    REQUIRE(x.size() == 5);
    REQUIRE(__D_EQ9(x[0].real(), 1.0));
    REQUIRE(__D_EQ9(x[0].imag(), 0.0));
    REQUIRE(__D_EQ9(x[1].real(), 2.0));
    REQUIRE(__D_EQ9(x[1].imag(), 0.0));
    REQUIRE(__D_EQ9(x[2].real(), 3.0));
    REQUIRE(__D_EQ9(x[2].imag(), 0.0));
    REQUIRE(__D_EQ9(x[3].real(), 4.0));
    REQUIRE(__D_EQ9(x[3].imag(), 0.0));
    REQUIRE(__D_EQ9(x[4].real(), 5.0));
    REQUIRE(__D_EQ9(x[4].imag(), 0.0));

    c.resize(16, 1);
    c << 32, -48, -8, 28, -8, 16, -16, 12, -16, 6, 10, -17, 10, 2, -4, 1;
    x = iexp::poly::complex_solve(c);
    REQUIRE(x.size() == 15);
    REQUIRE(__D_EQ7(x[0].real(), -1.6078107423472359));
    REQUIRE(__D_EQ7(x[0].imag(), 0.0));
    REQUIRE(__D_EQ7(x[1].real(), -1.3066982484920768));
    REQUIRE(__D_EQ7(x[1].imag(), 0.0));
    REQUIRE(__D_EQ7(x[2].real(), -1.0000000000000000));
    REQUIRE(__D_EQ7(x[2].imag(), 0.0));
    REQUIRE(__D_EQ7(x[4].real(), -0.65893856175240950));
    REQUIRE(__D_EQ7(x[4].imag(), -0.83459757287426684));
    REQUIRE(__D_EQ7(x[3].real(), -0.65893856175240950));
    REQUIRE(__D_EQ7(x[3].imag(), 0.83459757287426684));
    REQUIRE(__D_EQ7(x[14].real(), 2.0));
    REQUIRE(__D_EQ7(x[14].imag(), 0.0));

    // compile
    x = iexp::poly::complex_solve((c + c) * 2);
}
