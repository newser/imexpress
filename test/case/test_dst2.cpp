#include <catch.hpp>
#include <fft/dst2.h>
#include <fft/idst2.h>
#include <iostream>
#include <test_util.h>

using namespace Eigen;
using namespace std;

TEST_CASE("dst2_float")
{
    // r2r fwd
    Matrix<float, 3, 2, RowMajor> i_rr(3, 2), o_rr(3, 2), o_rr2(3, 2);

    i_rr << 1, 3, 2, 9, 4, 8;
    o_rr.array() = fft::dst2(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 77.297, 0.001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -38.9456, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), 1.08672, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 9.55177, 0.00001));

    // again
    o_rr.array() = fft::dst2(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 77.297, 0.001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -38.9456, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), 1.08672, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 9.55177, 0.00001));

    // inv
    o_rr2.array() = fft::idst2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 48, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 48, 0.0001));

    // again
    o_rr2.array() = fft::idst2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 48, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 48, 0.0001));

    // inv, normalize
    o_rr2.array() = fft::idst2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));

    // again
    o_rr2.array() = fft::idst2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));
}

TEST_CASE("dst2_double")
{
    // r2r fwd
    Matrix<double, 3, 2, RowMajor> i_rr(3, 2), o_rr(3, 2), o_rr2(3, 2);

    i_rr << 1, 3, 2, 9, 4, 8;
    o_rr.array() = fft::dst2(i_rr.array());
#if 0
    // octave's dst:
    o_rr *= sqrt(2.0/3) * sqrt(2.0/2)/4;
    o_rr.col(0) /= sqrt(2.0);
    o_rr.row(0) /= sqrt(2.0);
    cout << o_rr << endl;
#endif
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 77.297, 0.001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -38.9456, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), 1.08672, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 9.55177, 0.00001));

    // again
    o_rr.array() = fft::dst2(i_rr.array());
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 77.297, 0.001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -38.9456, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(2, 0), 1.08672, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr(2, 1), 9.55177, 0.00001));

    // inv
    o_rr2.array() = fft::idst2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 48, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 48, 0.0001));

    // again
    o_rr2.array() = fft::idst2(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3 * 48, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4 * 48, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8 * 48, 0.0001));

    // inv, normalize
    o_rr2.array() = fft::idst2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));

    // again
    o_rr2.array() = fft::idst2<true>(o_rr.array());
    REQUIRE(__F_EQ_IN(o_rr2(0, 0), 1, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(0, 1), 3, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 0), 4, 0.00001));
    REQUIRE(__F_EQ_IN(o_rr2(2, 1), 8, 0.00001));
}
