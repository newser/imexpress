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
    Matrix<float, Dynamic, Dynamic, RowMajor> i_r(2, 5);
    Matrix<std::complex<float>, Dynamic, Dynamic, RowMajor> o_r2(2, 5),
        o_r3(2, 5);

    i_r << 1, 9, 8, 7, 0, 9, 3, 1, 5, 4;
    o_r2.array() = fft::fft2(i_r.array());
    REQUIRE(__F_EQ8(o_r2(0, 0).real(), 47));
    REQUIRE(__F_EQ8(o_r2(0, 0).imag(), 0));
    REQUIRE(__F_EQ8(o_r2(0, 1).real(), -2.04508495));
    REQUIRE(__F_EQ8(o_r2(0, 1).imag(), -5.84509754));
    REQUIRE(__F_EQ8(o_r2(0, 4).real(), -2.04508495));
    REQUIRE(__F_EQ8(o_r2(0, 4).imag(), 5.84509754));
    REQUIRE(__F_EQ8(o_r2(1, 0).real(), 3));
    REQUIRE(__F_EQ8(o_r2(1, 0).imag(), 0));
    REQUIRE(__F_EQ7(o_r2(1, 4).real(), -14.6631193));
    REQUIRE(__F_EQ7(o_r2(1, 4).imag(), 12.4494925));

    // again
    o_r3 = fft::fft2(i_r);
    REQUIRE(__F_EQ8(o_r3(0, 0).real(), 47));
    REQUIRE(__F_EQ8(o_r3(0, 0).imag(), 0));
    REQUIRE(__F_EQ8(o_r3(0, 1).real(), -2.04508495));
    REQUIRE(__F_EQ8(o_r3(0, 1).imag(), -5.84509754));
    REQUIRE(__F_EQ8(o_r3(0, 4).real(), -2.04508495));
    REQUIRE(__F_EQ8(o_r3(0, 4).imag(), 5.84509754));
    REQUIRE(__F_EQ8(o_r3(1, 0).real(), 3));
    REQUIRE(__F_EQ8(o_r3(1, 0).imag(), 0));
    REQUIRE(__F_EQ7(o_r3(1, 4).real(), -14.6631193));
    REQUIRE(__F_EQ7(o_r3(1, 4).imag(), 12.4494925));
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

    Matrix<std::complex<double>, 4, 3> c;
    Matrix<std::complex<double>, 4, 3> c_o;
    c << 1, 2, 3, 4, 5, 6, 7, 7, 7, 9, 8, 7;
    c_o = fft::fft2(c);
    REQUIRE(__F_EQ7(c_o(0, 0).real(), 66));
    REQUIRE(__F_EQ7(c_o(0, 0).imag(), 0));
    REQUIRE(__F_EQ7(c_o(3, 0).real(), -15));
    REQUIRE(__F_EQ7(c_o(3, 0).imag(), -9));
    REQUIRE(__F_EQ7(c_o(3, 1).real(), -3.2320508075688772));
    REQUIRE(__F_EQ7(c_o(3, 1).imag(), -2.1339745962155616));

    // again
    c_o = fft::fft2(c);
    REQUIRE(__F_EQ7(c_o(0, 0).real(), 66));
    REQUIRE(__F_EQ7(c_o(0, 0).imag(), 0));
    REQUIRE(__F_EQ7(c_o(3, 0).real(), -15));
    REQUIRE(__F_EQ7(c_o(3, 0).imag(), -9));
    REQUIRE(__F_EQ7(c_o(3, 1).real(), -3.2320508075688772));
    REQUIRE(__F_EQ7(c_o(3, 1).imag(), -2.1339745962155616));

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

    // r2c
    {
        Matrix<double, 4, 3> c;
        Matrix<std::complex<double>, 4, 3> c_o;
        c << 1, 2, 3, 4, 5, 6, 7, 7, 7, 9, 8, 7;
        c_o = fft::fft2(c);
        REQUIRE(__F_EQ7(c_o(0, 0).real(), 66));
        REQUIRE(__F_EQ7(c_o(0, 0).imag(), 0));
        REQUIRE(__F_EQ7(c_o(3, 0).real(), -15));
        REQUIRE(__F_EQ7(c_o(3, 0).imag(), -9));
        REQUIRE(__F_EQ7(c_o(3, 1).real(), -3.2320508075688772));
        REQUIRE(__F_EQ7(c_o(3, 1).imag(), -2.1339745962155616));

        // again
        c_o = fft::fft2(c);
        REQUIRE(__F_EQ7(c_o(0, 0).real(), 66));
        REQUIRE(__F_EQ7(c_o(0, 0).imag(), 0));
        REQUIRE(__F_EQ7(c_o(3, 0).real(), -15));
        REQUIRE(__F_EQ7(c_o(3, 0).imag(), -9));
        REQUIRE(__F_EQ7(c_o(3, 1).real(), -3.2320508075688772));
        REQUIRE(__F_EQ7(c_o(3, 1).imag(), -2.1339745962155616));
    }
}
