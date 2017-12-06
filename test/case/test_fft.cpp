#include <catch.hpp>
#include <fft/fft.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("fft_float")
{
    VectorXcf i(8), o(8), o2(8);

    // c2c fwd
    i << 0, complex<float>(1, 1), complex<float>(2, 2), complex<float>(3, 3),
        complex<float>(4, 4), complex<float>(5, 5), complex<float>(6, 6),
        complex<float>(7, 7);
    o.array() = fft::fft(i.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -13.6569, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].real(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -13.6569, 0.0001));

    // again
    o.array() = fft::fft(i.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -13.6569, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].real(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -13.6569, 0.0001));

    // r2c
    VectorXf i_r(8), o_r2(8);

    i_r << 0, 1, 2, 3, 4, 5, 6, 7;
    o.array() = fft::fft(i_r.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -4, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 9.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[7].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -9.65685, 0.0001));

    // again
    o.array() = fft::fft(i_r.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -4, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 9.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[7].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -9.65685, 0.0001));
}

TEST_CASE("fft_double")
{
    VectorXcd i(8), o(8), o2(8);

    // c2c fwd
    i << 0, complex<double>(1, 1), complex<double>(2, 2), complex<double>(3, 3),
        complex<double>(4, 4), complex<double>(5, 5), complex<double>(6, 6),
        complex<double>(7, 7);
    o.array() = fft::fft(i.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -13.6569, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].real(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -13.6569, 0.0001));

    // again
    o.array() = fft::fft(i.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -13.6569, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].real(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -13.6569, 0.0001));

    // r2c
    VectorXd i_r(8), o_r2(8);

    i_r << 0, 1, 2, 3, 4, 5, 6, 7;
    o.array() = fft::fft(i_r.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -4, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 9.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[7].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -9.65685, 0.0001));

    // again
    o.array() = fft::fft(i_r.array());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -4, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 9.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[7].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -9.65685, 0.0001));
}
