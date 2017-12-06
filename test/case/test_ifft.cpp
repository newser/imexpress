#include <catch.hpp>
#include <fft/fft.h>
#include <fft/ifft.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("ifft_float")
{
    VectorXcf i(8), o(8), o2(8);

    // c2c fwd
    i << 0, complex<float>(1, 1), complex<float>(2, 2), complex<float>(3, 3),
        complex<float>(4, 4), complex<float>(5, 5), complex<float>(6, 6),
        complex<float>(7, 7);
    o.array() = fft::fft(i.array());

    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7 * 8, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7, 0.0001));

    // again
    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7 * 8, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7, 0.0001));

    // r2c
    VectorXf i_r(8), o_r2(8);

    i_r << 0, 1, 2, 3, 4, 5, 6, 7;
    o.array() = fft::fft(i_r.array());
    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));

    // again
    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));
}

TEST_CASE("ifft_double")
{
    VectorXcd i(8), o(8), o2(8);

    // c2c fwd
    i << 0, complex<double>(1, 1), complex<double>(2, 2), complex<double>(3, 3),
        complex<double>(4, 4), complex<double>(5, 5), complex<double>(6, 6),
        complex<double>(7, 7);
    o.array() = fft::fft(i.array());

    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7 * 8, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7, 0.0001));

    // again
    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7 * 8, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7, 0.0001));

    // r2c
    VectorXd i_r(8), o_r2(8);

    i_r << 0, 1, 2, 3, 4, 5, 6, 7;
    o.array() = fft::fft(i_r.array());
    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));

    // again
    o2.array() = fft::ifft(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));

    o2.array() = fft::ifft<true>(o.array());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7, 0.0001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 0, 0.0001));
}
