#include <catch.hpp>
#include <fft/dst.h>
#include <fft/idst.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("dst_float")
{
    // r2r fwd
    VectorXf i_rr(8), o_rr(8), o_rr2(8);

    i_rr << 0, 1, 2, 3, 4, 5, 6, 7;
    o_rr.array() = fft::dst(i_rr.array());
#if 0
    // octave's dst:
    o_rr /= 2;
#endif
    REQUIRE(__F_EQ_IN(o_rr[0], 39.699, 0.001));
    REQUIRE(__F_EQ_IN(o_rr[1], -24.7273, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -1.58694, 0.00001));

    // again
    o_rr.array() = fft::dst(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr[0], 39.699, 0.001));
    REQUIRE(__F_EQ_IN(o_rr[1], -24.7273, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -1.58694, 0.00001));

    // inv
    o_rr2.array() = fft::idst(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1 * 18, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7 * 18, 0.00001));

    // again
    o_rr2.array() = fft::idst(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1 * 18, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7 * 18, 0.00001));

    // inv, normalize
    o_rr2.array() = fft::idst<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7, 0.00001));

    // again
    o_rr2.array() = fft::idst<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7, 0.00001));
}

TEST_CASE("dst_double")
{
    // r2r fwd
    VectorXd i_rr(8), o_rr(8), o_rr2(8);

    i_rr << 0, 1, 2, 3, 4, 5, 6, 7;
    o_rr.array() = fft::dst(i_rr.array());
#if 0
    // octave's dst:
    o_rr /= 2;
#endif
    REQUIRE(__F_EQ_IN(o_rr[0], 39.699, 0.001));
    REQUIRE(__F_EQ_IN(o_rr[1], -24.7273, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -1.58694, 0.00001));

    Matrix<double, 4, 2, RowMajor> im, om, im2;
    im << 0, 1, 2, 3, 4, 5, 6, 7;
    om = fft::dst(im);
    REQUIRE(__F_EQ_IN(om(0, 0), 39.699, 0.001));
    REQUIRE(__F_EQ_IN(om(0, 1), -24.7273, 0.0001));
    REQUIRE(__F_EQ_IN(om(3, 1), -1.58694, 0.00001));

    // again
    o_rr.array() = fft::dst(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr[0], 39.699, 0.001));
    REQUIRE(__F_EQ_IN(o_rr[1], -24.7273, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -1.58694, 0.00001));

    // inv
    o_rr2.array() = fft::idst(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1 * 18, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7 * 18, 0.00001));

    im2 = fft::dst(om);
    REQUIRE(__F_EQ_IN(im2(0, 0), 0, 0.001));
    REQUIRE(__F_EQ_IN(im2(0, 1), 1 * 18, 0.0001));
    REQUIRE(__F_EQ_IN(im2(3, 1), 7 * 18, 0.00001));

    // again
    o_rr2.array() = fft::idst(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1 * 18, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7 * 18, 0.00001));

    // inv, normalize
    o_rr2.array() = fft::idst<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7, 0.00001));

    // again
    o_rr2.array() = fft::idst<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7, 0.00001));
}
