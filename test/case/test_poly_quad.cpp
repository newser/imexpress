#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <poly/quadratic.h>

TEST_CASE("poly_quadratic")
{
    double x0 = 0, x1 = 0;
    auto s4 = iexp::poly::solve_quad(4.0, -20.0, 26.0);
    REQUIRE(std::isnan(std::get<0>(s4)));
    REQUIRE(std::isnan(std::get<1>(s4)));

    auto s3 = iexp::poly::solve_quad(1, 0, 0);
    REQUIRE(__D_EQ9(std::get<0>(s3), 0));
    REQUIRE(__D_EQ9(std::get<1>(s3), 0));

    auto s5 = iexp::poly::solve_quad(0, 1, 1);
    REQUIRE(__D_EQ9(std::get<0>(s5), -1));
    REQUIRE(std::isnan(std::get<1>(s5)));

    std::complex<double> z0, z1;
    auto zz = iexp::poly::complex_solve_quad(4.0, -20.0, 26.0);
    REQUIRE(__D_EQ_IN(std::get<0>(zz).real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<0>(zz).imag(), -0.5, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<1>(zz).real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<1>(zz).imag(), 0.5, 1e-9));

    auto xx = iexp::poly::solve_quad(4.0, -20.0, 25.0);
    REQUIRE(__D_EQ_IN(std::get<0>(xx), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<1>(xx), 2.5, 1e-9));

    iexp::poly::complex_solve_quad(4.0, -20.0, 25.0, z0, z1);
    REQUIRE(__D_EQ_IN(z0.real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(z0.imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(z1.real(), 2.5, 1e-9));
    REQUIRE(__D_EQ_IN(z1.imag(), 0, 1e-9));

    auto s1 = iexp::poly::solve_quad(4.0, -20.0, 21.0);
    REQUIRE(__D_EQ_IN(std::get<0>(s1), 1.5, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<1>(s1), 3.5, 1e-9));

    auto s2 = iexp::poly::complex_solve_quad(4.0, -20.0, 21.0);
    REQUIRE(__D_EQ_IN(std::get<0>(s2).real(), 1.5, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<0>(s2).imag(), 0, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<1>(s2).real(), 3.5, 1e-9));
    REQUIRE(__D_EQ_IN(std::get<1>(s2).imag(), 0, 1e-9));

    auto t1 = iexp::poly::complex_solve_quad(0, 0, 21.0);
    REQUIRE(std::isnan(std::get<0>(t1).real()));
    REQUIRE(std::isnan(std::get<0>(t1).imag()));
    REQUIRE(std::isnan(std::get<1>(t1).real()));
    REQUIRE(std::isnan(std::get<1>(t1).imag()));

    auto t2 = iexp::poly::complex_solve_quad(0, 1, 21.0);
    REQUIRE(__D_EQ9(std::get<0>(t2).real(), -21.0));
    REQUIRE(__D_EQ9(std::get<0>(t2).imag(), 0));
    REQUIRE(std::isnan(std::get<1>(t2).real()));
    REQUIRE(std::isnan(std::get<1>(t2).imag()));
}
