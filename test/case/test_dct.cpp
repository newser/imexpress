#include <catch.hpp>
#include <fft/dct.h>
#include <fft/idct.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("dct_float")
{
    // r2r fwd
    VectorXf i_rr(8), o_rr(8), o_rr2(8);

    i_rr << 0, 1, 2, 3, 4, 5, 6, 7;
    o_rr.array() = fft::dct(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr[0], 56., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[1], -25.7693, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -0.20281, 0.00001));
    // octave's dct:
    // cout << o_rr[0] * std::sqrt(1.0/8)/2<< endl;
    // cout << o_rr[1] * std::sqrt(2.0/8)/2<< endl;
    // cout << o_rr[7] * std::sqrt(2.0/8)/2<< endl;

    MatrixXf im(2, 4), om(2, 4), im2(2, 4);
    im(0, 0) = 0;
    im(1, 0) = 1;
    im(0, 1) = 2;
    im(1, 1) = 3;
    im(0, 2) = 4;
    im(1, 2) = 5;
    im(0, 3) = 6;
    im(1, 3) = 7;
    om = fft::dct(im);
    REQUIRE(__F_EQ_IN(om(0, 0), 56., 0.0001));
    REQUIRE(__F_EQ_IN(om(1, 0), -25.7693, 0.0001));
    REQUIRE(__F_EQ_IN(om(1, 3), -0.20281, 0.00001));

    // again
    o_rr.array() = fft::dct(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr[0], 56., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[1], -25.7693, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -0.20281, 0.00001));

    // inv
    o_rr2.array() = fft::idct(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1. * 16, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7. * 16, 0.00001));

    im2 = fft::idct(om);
    REQUIRE(__F_EQ_IN(im2(0, 0), 0., 0.0001));
    REQUIRE(__F_EQ_IN(im2(1, 0), 1. * 16, 0.0001));
    REQUIRE(__F_EQ_IN(im2(1, 3), 7. * 16, 0.00001));

    // again
    o_rr2.array() = fft::idct(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1. * 16, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7. * 16, 0.00001));

    // inv, normalize
    o_rr2.array() = fft::idct<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7., 0.00001));

    // again
    o_rr2.array() = fft::idct<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7., 0.00001));
}

TEST_CASE("dct_double")
{
    // r2r fwd
    VectorXd i_rr(8), o_rr(8), o_rr2(8);

    i_rr << 0, 1, 2, 3, 4, 5, 6, 7;
    o_rr.array() = fft::dct(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr[0], 56., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[1], -25.7693, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -0.20281, 0.00001));
    // octave's dct:
    // cout << o_rr[0] * std::sqrt(1.0/8)/2<< endl;
    // cout << o_rr[1] * std::sqrt(2.0/8)/2<< endl;
    // cout << o_rr[7] * std::sqrt(2.0/8)/2<< endl;

    // again
    o_rr.array() = fft::dct(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr[0], 56., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[1], -25.7693, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -0.20281, 0.00001));

    // inv
    o_rr2.array() = fft::idct(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1. * 16, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7. * 16, 0.00001));

    // again
    o_rr2.array() = fft::idct(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1. * 16, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7. * 16, 0.00001));

    // inv, normalize
    o_rr2.array() = fft::idct<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7., 0.00001));

    // again
    o_rr2.array() = fft::idct<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2[0], 0., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[1], 1., 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2[7], 7., 0.00001));
}
