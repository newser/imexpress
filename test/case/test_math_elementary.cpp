#include <../test/test_util.h>
#include <catch.hpp>
#include <gsl/gsl_math.h>
#include <iostream>
#include <math/elementary.h>
#include <unsupported/Eigen/SpecialFunctions>

TEST_CASE("math_elementary")
{
    SECTION("exp, log, log1p, log10")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random();
        iexp::Matrix2d m2 = iexp::log1p(m.array());
        REQUIRE(::gsl_log1p(m(0, 0)) == m2(0, 0));
        REQUIRE(__D_EQ9(::gsl_log1p(m(0, 1)), m2(0, 1)));
        REQUIRE(__D_EQ9(::gsl_log1p(m(1, 0)), m2(1, 0)));
        REQUIRE(::gsl_log1p(m(1, 1)) == m2(1, 1));

        m2 = iexp::exp(m.array());
        REQUIRE(std::exp(m(0, 0)) == m2(0, 0));
        REQUIRE(std::exp(m(0, 1)) == m2(0, 1));
        REQUIRE(std::exp(m(1, 0)) == m2(1, 0));
        REQUIRE(std::exp(m(1, 1)) == m2(1, 1));

        m2 = iexp::log(m.cwiseAbs().array());
        REQUIRE(std::log(std::abs(m(0, 0))) == m2(0, 0));
        REQUIRE(std::log(std::abs(m(0, 1))) == m2(0, 1));
        REQUIRE(std::log(std::abs(m(1, 0))) == m2(1, 0));
        REQUIRE(std::log(std::abs(m(1, 1))) == m2(1, 1));

        m2 = iexp::log10(m.cwiseAbs().array());
        REQUIRE(std::log10(std::abs(m(0, 0))) == m2(0, 0));
        REQUIRE(std::log10(std::abs(m(0, 1))) == m2(0, 1));
        REQUIRE(std::log10(std::abs(m(1, 0))) == m2(1, 0));
        REQUIRE(std::log10(std::abs(m(1, 1))) == m2(1, 1));
    }

    SECTION("real, imag, conj")
    {
        iexp::Matrix2cd m = iexp::Matrix2cd::Random();
        iexp::Matrix2cd r = iexp::real(m.array());
        iexp::Matrix2cd i = iexp::imag(m.array());
        iexp::Matrix2cd c = iexp::conj(m.array());
        REQUIRE(std::real(m(0, 0)) == r(0, 0));
        REQUIRE(std::real(m(0, 1)) == r(0, 1));
        REQUIRE(std::real(m(1, 0)) == r(1, 0));
        REQUIRE(std::real(m(1, 1)) == r(1, 1));
        REQUIRE(std::imag(m(0, 0)) == i(0, 0));
        REQUIRE(std::imag(m(0, 1)) == i(0, 1));
        REQUIRE(std::imag(m(1, 0)) == i(1, 0));
        REQUIRE(std::imag(m(1, 1)) == i(1, 1));
        REQUIRE(std::conj(m(0, 0)) == c(0, 0));
        REQUIRE(std::conj(m(0, 1)) == c(0, 1));
        REQUIRE(std::conj(m(1, 0)) == c(1, 0));
        REQUIRE(std::conj(m(1, 1)) == c(1, 1));
    }

    SECTION("abs, inverse, conj")
    {
        iexp::Matrix2cd m = iexp::Matrix2cd::Random();
        iexp::Matrix2cd a = iexp::abs(m.array());
        iexp::Matrix2cd i = iexp::inverse(m.array());
        REQUIRE(std::abs(m(0, 0)) == a(0, 0));
        REQUIRE(std::abs(m(0, 1)) == a(0, 1));
        REQUIRE(std::abs(m(1, 0)) == a(1, 0));
        REQUIRE(std::abs(m(1, 1)) == a(1, 1));
    }

    SECTION("pow, sqrt, rsqrt, square, cube")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random();
        iexp::Matrix2d m2 = iexp::pow(m.array(), m.array());
        REQUIRE(__D_EQ(std::pow(m(0, 0), m(0, 0)), m2(0, 0)));
        REQUIRE(__D_EQ(std::pow(m(0, 1), m(0, 1)), m2(0, 1)));
        REQUIRE(__D_EQ(std::pow(m(1, 0), m(1, 0)), m2(1, 0)));
        REQUIRE(__D_EQ(std::pow(m(1, 1), m(1, 1)), m2(1, 1)));

        m2 = iexp::sqrt(m.array().abs());
        REQUIRE(__D_EQ(std::sqrt(std::abs(m(0, 0))), m2(0, 0)));
        REQUIRE(__D_EQ(std::sqrt(std::abs(m(0, 1))), m2(0, 1)));
        REQUIRE(__D_EQ(std::sqrt(std::abs(m(1, 0))), m2(1, 0)));
        REQUIRE(__D_EQ(std::sqrt(std::abs(m(1, 1))), m2(1, 1)));

        m2 = iexp::rsqrt(m.array().abs());
        REQUIRE(__D_EQ(1 / std::sqrt(std::abs(m(0, 0))), m2(0, 0)));
        REQUIRE(__D_EQ(1 / std::sqrt(std::abs(m(0, 1))), m2(0, 1)));
        REQUIRE(__D_EQ(1 / std::sqrt(std::abs(m(1, 0))), m2(1, 0)));
        REQUIRE(__D_EQ(1 / std::sqrt(std::abs(m(1, 1))), m2(1, 1)));

        m2 = iexp::square(m.array());
        REQUIRE(__D_EQ(m(0, 0) * m(0, 0), m2(0, 0)));
        REQUIRE(__D_EQ(m(0, 1) * m(0, 1), m2(0, 1)));
        REQUIRE(__D_EQ(m(1, 0) * m(1, 0), m2(1, 0)));
        REQUIRE(__D_EQ(m(1, 1) * m(1, 1), m2(1, 1)));

        m2 = iexp::cube(m.array());
        REQUIRE(__D_EQ(m(0, 0) * m(0, 0) * m(0, 0), m2(0, 0)));
        REQUIRE(__D_EQ(m(0, 1) * m(0, 1) * m(0, 1), m2(0, 1)));
        REQUIRE(__D_EQ(m(1, 0) * m(1, 0) * m(1, 0), m2(1, 0)));
        REQUIRE(__D_EQ(m(1, 1) * m(1, 1) * m(1, 1), m2(1, 1)));

        m2 = iexp::abs2(m.array());
        REQUIRE(__D_EQ(m(0, 0) * m(0, 0), m2(0, 0)));
        REQUIRE(__D_EQ(m(0, 1) * m(0, 1), m2(0, 1)));
        REQUIRE(__D_EQ(m(1, 0) * m(1, 0), m2(1, 0)));
        REQUIRE(__D_EQ(m(1, 1) * m(1, 1), m2(1, 1)));
    }

    SECTION("Trigonometric")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random();
        iexp::Matrix2d m2 = iexp::sin(m.array());
        REQUIRE(std::sin(m(0, 0)) == m2(0, 0));
        REQUIRE(std::sin(m(0, 1)) == m2(0, 1));
        REQUIRE(std::sin(m(1, 0)) == m2(1, 0));
        REQUIRE(std::sin(m(1, 1)) == m2(1, 1));

        m2 = iexp::cos(m.array());
        REQUIRE(std::cos(m(0, 0)) == m2(0, 0));
        REQUIRE(std::cos(m(0, 1)) == m2(0, 1));
        REQUIRE(std::cos(m(1, 0)) == m2(1, 0));
        REQUIRE(std::cos(m(1, 1)) == m2(1, 1));

        m2 = iexp::tan(m.array());
        REQUIRE(std::tan(m(0, 0)) == m2(0, 0));
        REQUIRE(std::tan(m(0, 1)) == m2(0, 1));
        REQUIRE(std::tan(m(1, 0)) == m2(1, 0));
        REQUIRE(std::tan(m(1, 1)) == m2(1, 1));

        m2 = iexp::asin(m.array());
        REQUIRE(std::asin(m(0, 0)) == m2(0, 0));
        REQUIRE(std::asin(m(0, 1)) == m2(0, 1));
        REQUIRE(std::asin(m(1, 0)) == m2(1, 0));
        REQUIRE(std::asin(m(1, 1)) == m2(1, 1));

        m2 = iexp::acos(m.array());
        REQUIRE(std::acos(m(0, 0)) == m2(0, 0));
        REQUIRE(std::acos(m(0, 1)) == m2(0, 1));
        REQUIRE(std::acos(m(1, 0)) == m2(1, 0));
        REQUIRE(std::acos(m(1, 1)) == m2(1, 1));

        m2 = iexp::atan(m.array());
        REQUIRE(std::atan(m(0, 0)) == m2(0, 0));
        REQUIRE(std::atan(m(0, 1)) == m2(0, 1));
        REQUIRE(std::atan(m(1, 0)) == m2(1, 0));
        REQUIRE(std::atan(m(1, 1)) == m2(1, 1));
    }

    SECTION("Hyperbolic")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random();
        iexp::Matrix2d m2 = iexp::sinh(m.array());
        REQUIRE(std::sinh(m(0, 0)) == m2(0, 0));
        REQUIRE(std::sinh(m(0, 1)) == m2(0, 1));
        REQUIRE(std::sinh(m(1, 0)) == m2(1, 0));
        REQUIRE(std::sinh(m(1, 1)) == m2(1, 1));

        m2 = iexp::cosh(m.array());
        REQUIRE(std::cosh(m(0, 0)) == m2(0, 0));
        REQUIRE(std::cosh(m(0, 1)) == m2(0, 1));
        REQUIRE(std::cosh(m(1, 0)) == m2(1, 0));
        REQUIRE(std::cosh(m(1, 1)) == m2(1, 1));

        m2 = iexp::tanh(m.array());
        REQUIRE(std::tanh(m(0, 0)) == m2(0, 0));
        REQUIRE(std::tanh(m(0, 1)) == m2(0, 1));
        REQUIRE(std::tanh(m(1, 0)) == m2(1, 0));
        REQUIRE(std::tanh(m(1, 1)) == m2(1, 1));
    }

    SECTION("Nearest integer")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random();
        iexp::Matrix2d m2 = iexp::ceil(m.array());
        REQUIRE(std::ceil(m(0, 0)) == m2(0, 0));
        REQUIRE(std::ceil(m(0, 1)) == m2(0, 1));
        REQUIRE(std::ceil(m(1, 0)) == m2(1, 0));
        REQUIRE(std::ceil(m(1, 1)) == m2(1, 1));

        m2 = iexp::floor(m.array());
        REQUIRE(std::floor(m(0, 0)) == m2(0, 0));
        REQUIRE(std::floor(m(0, 1)) == m2(0, 1));
        REQUIRE(std::floor(m(1, 0)) == m2(1, 0));
        REQUIRE(std::floor(m(1, 1)) == m2(1, 1));

        m2 = iexp::round(m.array());
        REQUIRE(std::round(m(0, 0)) == m2(0, 0));
        REQUIRE(std::round(m(0, 1)) == m2(0, 1));
        REQUIRE(std::round(m(1, 0)) == m2(1, 0));
        REQUIRE(std::round(m(1, 1)) == m2(1, 1));
    }

    SECTION("Nearest integer")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random();
        iexp::Matrix2d m2 = iexp::ceil(m.array());
        REQUIRE(std::ceil(m(0, 0)) == m2(0, 0));
        REQUIRE(std::ceil(m(0, 1)) == m2(0, 1));
        REQUIRE(std::ceil(m(1, 0)) == m2(1, 0));
        REQUIRE(std::ceil(m(1, 1)) == m2(1, 1));

        m2 = iexp::floor(m.array());
        REQUIRE(std::floor(m(0, 0)) == m2(0, 0));
        REQUIRE(std::floor(m(0, 1)) == m2(0, 1));
        REQUIRE(std::floor(m(1, 0)) == m2(1, 0));
        REQUIRE(std::floor(m(1, 1)) == m2(1, 1));

        m2 = iexp::round(m.array());
        REQUIRE(std::round(m(0, 0)) == m2(0, 0));
        REQUIRE(std::round(m(0, 1)) == m2(0, 1));
        REQUIRE(std::round(m(1, 0)) == m2(1, 0));
        REQUIRE(std::round(m(1, 1)) == m2(1, 1));
    }

    SECTION("Error and gamma functions")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random().cwiseAbs();
        iexp::Matrix2d m2 = iexp::erf(m.array());
        REQUIRE(std::erf(m(0, 0)) == m2(0, 0));
        REQUIRE(std::erf(m(0, 1)) == m2(0, 1));
        REQUIRE(std::erf(m(1, 0)) == m2(1, 0));
        REQUIRE(std::erf(m(1, 1)) == m2(1, 1));

        m2 = iexp::erfc(m.array());
        REQUIRE(std::erfc(m(0, 0)) == m2(0, 0));
        REQUIRE(std::erfc(m(0, 1)) == m2(0, 1));
        REQUIRE(std::erfc(m(1, 0)) == m2(1, 0));
        REQUIRE(std::erfc(m(1, 1)) == m2(1, 1));

        m2 = iexp::lgamma(m.array());
        REQUIRE(std::lgamma(m(0, 0)) == m2(0, 0));
        REQUIRE(std::lgamma(m(0, 1)) == m2(0, 1));
        REQUIRE(std::lgamma(m(1, 0)) == m2(1, 0));
        REQUIRE(std::lgamma(m(1, 1)) == m2(1, 1));
        m2 = iexp::lgamma(m.array());

        m2 = iexp::digamma(m.array());
        m2 = iexp::igamma(m.array(), m.array());
        m2 = iexp::igammac(m.array(), m.array());
        m2 = iexp::polygamma(m.array(), m.array());
        m2 = iexp::betainc(m.array(), m.array(), m.array());
        m2 = iexp::zeta(m.array(), m.array());
    }

    SECTION("Sign")
    {
        iexp::Matrix2d m = iexp::Matrix2d::Random();
        iexp::Matrix2d m2 = iexp::sign(m.array());
        std::cout << m2 << std::endl;
        REQUIRE(m(0, 0) * m2(0, 0) >= 0);
        REQUIRE(m(0, 1) * m2(0, 1) >= 0);
        REQUIRE(m(1, 0) * m2(1, 0) >= 0);
        REQUIRE(m(1, 1) * m2(1, 1) >= 0);
    }
}
