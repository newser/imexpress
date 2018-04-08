#include <catch.hpp>
#include <fft/fft2.h>
#include <fft/ifft2.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("ifft2_float")
{
    Matrix<complex<float>, Dynamic, Dynamic, RowMajor> i(3, 2), o(3, 2),
        o2(3, 2);

    // c2c inv
    i << complex<float>(1, 2), complex<float>(3, 4), complex<float>(8, 0),
        complex<float>(2, 2), complex<float>(9, 2), complex<float>(3, 9);
    o.array() = fft::fft2(i.array());

    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 2 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 3 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 9 * 6));

    // again
    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 2 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 3 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 9 * 6));

    // again, matrix col major
    Matrix<complex<float>, 3, 2> c = o;
    Matrix<complex<float>, 3, 2> c_o = fft::ifft2(c);
    REQUIRE(__F_EQ5(c_o(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(c_o(0, 0).imag(), 2 * 6));
    REQUIRE(__F_EQ5(c_o(2, 1).real(), 3 * 6));
    REQUIRE(__F_EQ5(c_o(2, 1).imag(), 9 * 6));

    // normalize
    o2.array() = fft::ifft2<true>(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 2));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 3));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 9));

    // r2c inv
    Matrix<float, Dynamic, Dynamic, RowMajor> i_r(3, 2), o_r2(3, 2);

    i_r << 1, 9, 8, 7, 0, 9;
    o.array() = fft::fft2(i_r.array());

    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 0 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 9 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 0 * 6));

    // again
    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 0 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 9 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 0 * 6));

    o2.array() = fft::ifft2<true>(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 0));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 9));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 0));
}

TEST_CASE("ifft2_double")
{
    Matrix<complex<double>, Dynamic, Dynamic, RowMajor> i(3, 2), o(3, 2),
        o2(3, 2);

    // c2c inv
    i << complex<double>(1, 2), complex<double>(3, 4), complex<double>(8, 0),
        complex<double>(2, 2), complex<double>(9, 2), complex<double>(3, 9);
    o.array() = fft::fft2(i.array());

    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 2 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 3 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 9 * 6));

    // again
    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 2 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 3 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 9 * 6));

    // normalize
    o2.array() = fft::ifft2<true>(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 2));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 3));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 9));

    // r2c inv
    Matrix<double, Dynamic, Dynamic, RowMajor> i_r(3, 2), o_r2(3, 2);

    i_r << 1, 9, 8, 7, 0, 9;
    o.array() = fft::fft2(i_r.array());

    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 0 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 9 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 0 * 6));

    // again
    o2.array() = fft::ifft2(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 0 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 9 * 6));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 0 * 6));

    o2.array() = fft::ifft2<true>(o.array());
    REQUIRE(__F_EQ5(o2(0, 0).real(), 1));
    REQUIRE(__F_EQ5(o2(0, 0).imag(), 0));
    REQUIRE(__F_EQ5(o2(2, 1).real(), 9));
    REQUIRE(__F_EQ5(o2(2, 1).imag(), 0));

    // again, matrix col major
    Matrix<std::complex<double>, 3, 2> c = o;
    Matrix<std::complex<double>, 3, 2> c_o = fft::ifft2(c);
    REQUIRE(__F_EQ5(c_o(0, 0).real(), 1 * 6));
    REQUIRE(__F_EQ5(c_o(0, 0).imag(), 0));
    REQUIRE(__F_EQ5(c_o(2, 1).real(), 9 * 6));
    REQUIRE(__F_EQ5(c_o(2, 1).imag(), 0));
}
