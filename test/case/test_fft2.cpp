#include <catch.hpp>
#include <fft/fft2.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("fft2_float")
{
    Matrix<complex<float>, Dynamic, Dynamic, RowMajor> i(3, 3), o(3, 3),
        o2(3, 3);

    // c2c fwd
    i << complex<float>(1, 2), complex<float>(3, 4), complex<float>(8, 0),
        complex<float>(2, 2), complex<float>(9, 2), complex<float>(3, 9),
        complex<float>(1, 1), complex<float>(3, 6), complex<float>(7, 7);
    o.array() = fft::fft2(i.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 37));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 33));
    REQUIRE(__F_EQ5(o(1, 2).real(), 1.80385));
    REQUIRE(__F_EQ5(o(1, 2).imag(), -5.19615));
    REQUIRE(__F_EQ5(o(2, 2).real(), -16.66025));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 3.80385));

    // again
    o.array() = fft::fft2(i.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 37));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 33));
    REQUIRE(__F_EQ5(o(1, 2).real(), 1.80385));
    REQUIRE(__F_EQ5(o(1, 2).imag(), -5.19615));
    REQUIRE(__F_EQ5(o(2, 2).real(), -16.66025));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 3.80385));

    // r2c
    Matrix<float, Dynamic, Dynamic, RowMajor> i_r(3, 3), o_r2(3, 3);

    i_r << 1, 9, 8, 7, 0, 9, 3, 1, 5;
    o.array() = fft::fft2(i_r.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 43));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 0));
    REQUIRE(__F_EQ5(o(1, 2).real(), -12.5));
    REQUIRE(__F_EQ5(o(1, 2).imag(), 4.33013));
    REQUIRE(__F_EQ5(o(2, 2).real(), -5));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 8.66025));

    // again
    o.array() = fft::fft2(i_r.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 43));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 0));
    REQUIRE(__F_EQ5(o(1, 2).real(), -12.5));
    REQUIRE(__F_EQ5(o(1, 2).imag(), 4.33013));
    REQUIRE(__F_EQ5(o(2, 2).real(), -5));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 8.66025));
}

TEST_CASE("fft2_double")
{
    Matrix<complex<double>, Dynamic, Dynamic, RowMajor> i(3, 3), o(3, 3),
        o2(3, 3);

    // c2c fwd
    i << complex<double>(1, 2), complex<double>(3, 4), complex<double>(8, 0),
        complex<double>(2, 2), complex<double>(9, 2), complex<double>(3, 9),
        complex<double>(1, 1), complex<double>(3, 6), complex<double>(7, 7);
    o.array() = fft::fft2(i.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 37));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 33));
    REQUIRE(__F_EQ5(o(1, 2).real(), 1.80385));
    REQUIRE(__F_EQ5(o(1, 2).imag(), -5.19615));
    REQUIRE(__F_EQ5(o(2, 2).real(), -16.66025));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 3.80385));

    // again
    o.array() = fft::fft2(i.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 37));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 33));
    REQUIRE(__F_EQ5(o(1, 2).real(), 1.80385));
    REQUIRE(__F_EQ5(o(1, 2).imag(), -5.19615));
    REQUIRE(__F_EQ5(o(2, 2).real(), -16.66025));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 3.80385));

    // r2c
    Matrix<double, Dynamic, Dynamic, RowMajor> i_r(3, 3), o_r2(3, 3);

    i_r << 1, 9, 8, 7, 0, 9, 3, 1, 5;
    o.array() = fft::fft2(i_r.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 43));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 0));
    REQUIRE(__F_EQ5(o(1, 2).real(), -12.5));
    REQUIRE(__F_EQ5(o(1, 2).imag(), 4.33013));
    REQUIRE(__F_EQ5(o(2, 2).real(), -5));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 8.66025));

    // again
    o.array() = fft::fft2(i_r.array());
    REQUIRE(__F_EQ5(o(0, 0).real(), 43));
    REQUIRE(__F_EQ5(o(0, 0).imag(), 0));
    REQUIRE(__F_EQ5(o(1, 2).real(), -12.5));
    REQUIRE(__F_EQ5(o(1, 2).imag(), 4.33013));
    REQUIRE(__F_EQ5(o(2, 2).real(), -5));
    REQUIRE(__F_EQ5(o(2, 2).imag(), 8.66025));
}
