#include <../test/test_util.h>
#include <catch.hpp>
#include <complex/complex.h>
#include <iostream>
#include <math/constant.h>

TEST_CASE("complex")
{
    std::complex<double> a(1, 2);
    REQUIRE(__D_EQ(std::real(a), 1));
    REQUIRE(__D_EQ(std::imag(a), 2));

    auto b = std::polar<double>(1, IEXP_PI_4);
    REQUIRE(__D_EQ(std::real(b), std::cos(IEXP_PI_4)));
    REQUIRE(__D_EQ(std::imag(b), std::sin(IEXP_PI_4)));

    a.real(std::cos(IEXP_PI_4));
    a.imag(std::sin(IEXP_PI_4));
    REQUIRE(a == b);
    REQUIRE(__D_EQ(std::arg(a), IEXP_PI_4));
    REQUIRE(__D_EQ(std::abs(a), 1));
    REQUIRE(__D_EQ(std::norm(a), 1));

    std::complex<double> c = a + b;
    REQUIRE(__D_EQ(c.real(), a.real() + b.real()));
    REQUIRE(__D_EQ(c.imag(), a.imag() + b.imag()));

    c = a - b;
    REQUIRE(__D_EQ(c.real(), a.real() - b.real()));
    REQUIRE(__D_EQ(c.imag(), a.imag() - b.imag()));

    c = a * b;
    REQUIRE(__D_EQ(c.real(), a.real() * b.real() - a.imag() * b.imag()));
    REQUIRE(__D_EQ(c.imag(), a.imag() * b.real() + a.real() * b.imag()));

    c = a + 1.;
    REQUIRE(__D_EQ(c.real(), a.real() + 1.));

    c = a + std::complex<double>(0, 1);
    REQUIRE(__D_EQ(c.imag(), a.imag() + 1.));

    c = std::conj(a);
    REQUIRE(__D_EQ(c.imag(), -a.imag()));

    c = std::sqrt(a);
    REQUIRE(__D_EQ(std::arg(c), IEXP_PI_4 / 2));
}
