#include <catch.hpp>
#include <fft/dct2.h>
#include <fft/idct2.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("dct2_float")
{
    // r2r fwd
    Matrix<float, 3, 2, RowMajor> i_rr(3, 2), o_rr(3, 2), o_rr2(3, 2);

    i_rr << 1, 3, 2, 9, 4, 8;
    o_rr.array() = fft::dct2(i_rr.array());
#if 0
    // octave's dct:
    o_rr *= sqrt(2.0/3) * sqrt(2.0/2)/4;
    o_rr.col(0) /= sqrt(2.0);
    o_rr.row(0) /= sqrt(2.0);
    cout << o_rr << endl;
#endif
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 108, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -36.7696, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), -12, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 11.3137, 0.00001));

    // again
    o_rr.array() = fft::dct2(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 108, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -36.7696, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), -12, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 11.3137, 0.00001));

    // inv
    o_rr2.array() = fft::idct2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 24, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 24, 0.00001));

    // again
    o_rr2.array() = fft::idct2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 24, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 24, 0.00001));

    // inv, normalize
    o_rr2.array() = fft::idct2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));

    // again
    o_rr2.array() = fft::idct2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));

    Matrix<float, 3, 2> ic_rr = o_rr;
    Matrix<float, 3, 2> oc_rr = fft::idct2<true>(ic_rr);
    REQUIRE(__F_EQ_IN(oc_rr(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(oc_rr(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(oc_rr(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(oc_rr(2, 1), 8, 0.00001));
}

TEST_CASE("dct2_double")
{
    // r2r fwd
    Matrix<double, 3, 2, RowMajor> i_rr(3, 2), o_rr(3, 2), o_rr2(3, 2);

    i_rr << 1, 3, 2, 9, 4, 8;
    o_rr.array() = fft::dct2(i_rr.array());
#if 0
    // octave's dct:
    o_rr *= sqrt(2.0/3) * sqrt(2.0/2)/4;
    o_rr.col(0) /= sqrt(2.0);
    o_rr.row(0) /= sqrt(2.0);
    cout << o_rr << endl;
#endif
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 108, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -36.7696, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), -12, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 11.3137, 0.00001));

    // again
    o_rr.array() = fft::dct2(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 108, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -36.7696, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), -12, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 11.3137, 0.00001));

    Matrix<double, 3, 2> ic_rr = i_rr;
    Matrix<double, 3, 2> oc_rr = fft::dct2(ic_rr);
    REQUIRE(__F_EQ_IN(oc_rr(0, 0), 108, 0.0001));
    REQUIRE(__F_EQ_IN(oc_rr(0, 1), -36.7696, 0.0001));
    REQUIRE(__F_EQ_IN(oc_rr(2, 0), -12, 0.00001));
    REQUIRE(__F_EQ_IN(oc_rr(2, 1), 11.3137, 0.00001));

    // inv
    o_rr2.array() = fft::idct2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 24, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 24, 0.00001));

    // again
    o_rr2.array() = fft::idct2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 24, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 24, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 24, 0.00001));

    // inv, normalize
    o_rr2.array() = fft::idct2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));

    // again
    o_rr2.array() = fft::idct2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));
}
