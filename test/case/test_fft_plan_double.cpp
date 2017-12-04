#include <catch.hpp>
#include <fft/fftw/plan_double.h>
#include <iostream>
#include <math/constant.h>
#include <test_util.h>

using namespace iexp;
using namespace std;

TEST_CASE("test_fft_plan_double_1d")
{
    fftw3::plan<double> p, q;
    VectorXcd i(8), o(8), o2(8);

    // c2c fwd
    i << 0, complex<double>(1, 1), complex<double>(2, 2), complex<double>(3, 3),
        complex<double>(4, 4), complex<double>(5, 5), complex<double>(6, 6),
        complex<double>(7, 7);
    p.fwd(i.size(), i.data(), o.data());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -13.6569, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].real(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -13.6569, 0.0001));

    // reuse
    p.fwd(i.size(), i.data(), o.data());
    REQUIRE(__F_EQ_IN(o[7].real(), 5.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[7].imag(), -13.6569, 0.0001));

    // c2c inv
    q.inv(i.size(), o.data(), o2.data());
    REQUIRE(__F_EQ_IN(o2[0].real(), 0 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[0].imag(), 0 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[1].real(), 1 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o2[1].imag(), 1 * 8, 0.00001));
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.00001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7 * 8, 0.0001));

    q.inv(i.size(), o.data(), o2.data());
    REQUIRE(__F_EQ_IN(o2[7].real(), 7 * 8, 0.00001));
    REQUIRE(__F_EQ_IN(o2[7].imag(), 7 * 8, 0.0001));

    // r2c fwd
    fftw3::plan<double> p_r, q_r;
    VectorXd i_r(8), o_r2(8);

    i_r << 0, 1, 2, 3, 4, 5, 6, 7;
    p_r.fwd(i_r.size(), i_r.data(), o.data());
    REQUIRE(__F_EQ_IN(o[0].real(), 28, 0.0001));
    REQUIRE(__F_EQ_IN(o[0].imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].real(), -4, 0.0001));
    REQUIRE(__F_EQ_IN(o[1].imag(), 9.65685, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].imag(), 0, 0.0001));

    p_r.fwd(i_r.size(), i_r.data(), o.data());
    REQUIRE(__F_EQ_IN(o[4].real(), -4, 0.00001));
    REQUIRE(__F_EQ_IN(o[4].imag(), 0, 0.0001));

    // c2r inv
    q_r.inv(o.size(), o.data(), o_r2.data());
    REQUIRE(__F_EQ_IN(o_r2[0], 0 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o_r2[1], 1 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o_r2[7], 7 * 8, 0.00001));

    q_r.inv(o.size(), o.data(), o_r2.data());
    REQUIRE(__F_EQ_IN(o_r2[7], 7 * 8, 0.00001));

    // r2r fwd
    fftw3::plan<double> p_rr, q_rr;
    VectorXd i_rr(8), o_rr(8), o_rr2(8);

    i_rr << 0, 1, 2, 3, 4, 5, 6, 7;
    p_rr.fwd<DCT_II>(i_rr.size(), i_rr.data(), o_rr.data());
    REQUIRE(__F_EQ_IN(o_rr[0], 56, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[1], -25.7693, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr[7], -0.20281, 0.00001));
    // octave's dct:
    // cout << o_rr[0] * std::sqrt(1.0/8)/2<< endl;
    // cout << o_rr[1] * std::sqrt(2.0/8)/2<< endl;
    // cout << o_rr[7] * std::sqrt(2.0/8)/2<< endl;

    p_rr.fwd<DCT_II>(i_rr.size(), i_rr.data(), o_rr.data());
    REQUIRE(__F_EQ_IN(o_rr[7], -0.20281, 0.00001));

    // r2r inv
    q_rr.inv<DCT_II>(o_rr.size(), o_rr.data(), o_rr2.data());
    REQUIRE(__F_EQ_IN(o_r2[0], 0 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o_r2[1], 1 * 8, 0.0001));
    REQUIRE(__F_EQ_IN(o_r2[7], 7 * 8, 0.00001));

    q_rr.inv<DCT_II>(o_rr.size(), o_rr.data(), o_rr2.data());
    REQUIRE(__F_EQ_IN(o_r2[7], 7 * 8, 0.00001));
}

TEST_CASE("test_fft_plan_double_2d")
{
    // c2c 2d
    fftw3::plan<double> p, q;
    Matrix<complex<double>, 3, 3, RowMajor> i(3, 3), o(3, 3), o2(3, 3),
        save_o(3, 3);

    i << 0, complex<double>(0, 1), complex<double>(0, 2), complex<double>(1, 0),
        complex<double>(1, 1), complex<double>(1, 2), complex<double>(2, 0),
        complex<double>(2, 1), complex<double>(2, 2);
    p.fwd(3, 3, i.data(), o.data());
    REQUIRE(__F_EQ_IN(o(0, 0).real(), 9, 0.0001));
    REQUIRE(__F_EQ_IN(o(0, 0).imag(), 9, 0.0001));
    REQUIRE(__F_EQ_IN(o(2, 0).real(), -4.5, 0.0001));
    REQUIRE(__F_EQ_IN(o(2, 0).imag(), -2.598076, 0.0001));

    p.fwd(3, 3, i.data(), o.data());
    REQUIRE(__F_EQ_IN(o(0, 0).real(), 9, 0.0001));
    REQUIRE(__F_EQ_IN(o(0, 0).imag(), 9, 0.0001));
    REQUIRE(__F_EQ_IN(o(2, 0).real(), -4.5, 0.0001));
    REQUIRE(__F_EQ_IN(o(2, 0).imag(), -2.598076, 0.0001));

    q.inv(3, 3, o.data(), o2.data());
    REQUIRE(__F_EQ_IN(o2(0, 0).real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2(0, 0).imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2(2, 0).real(), 2 * 9, 0.0001));
    REQUIRE(__F_EQ_IN(o2(2, 0).imag(), 0, 0.0001));

    q.inv(3, 3, o.data(), o2.data());
    REQUIRE(__F_EQ_IN(o2(0, 0).real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2(0, 0).imag(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(o2(2, 0).real(), 2 * 9, 0.0001));
    REQUIRE(__F_EQ_IN(o2(2, 0).imag(), 0, 0.0001));

    // r2c 2d
    fftw3::plan<double> p2, q2;
    Matrix<double, 3, 3, RowMajor> i_r(3, 3), o2_r(3, 3);
    Matrix<complex<double>, 3, 2, RowMajor> oo_r(3, 2);

    i_r << 0, 1, 2, 3, 4, 5, 6, 7, 8;
    oo_r.fill(complex<double>(9, 9));
    p2.fwd(3, 3, i_r.data(), oo_r.data());
    REQUIRE(__F_EQ_IN(oo_r(0, 1).real(), -4.5, 0.0001));
    REQUIRE(__F_EQ_IN(oo_r(0, 1).imag(), 2.59808, 0.00001));
    REQUIRE(__F_EQ_IN(oo_r(2, 0).real(), -13.5, 0.0001));
    REQUIRE(__F_EQ_IN(oo_r(2, 0).imag(), -7.79423, 0.0001));
    REQUIRE(__F_EQ_IN(oo_r(2, 1).real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(oo_r(2, 1).imag(), 0, 0.0001));

    p2.fwd(3, 3, i_r.data(), o.data());
    REQUIRE(__F_EQ_IN(oo_r(2, 0).real(), -13.5, 0.0001));
    REQUIRE(__F_EQ_IN(oo_r(2, 0).imag(), -7.79423, 0.0001));
    REQUIRE(__F_EQ_IN(oo_r(2, 1).real(), 0, 0.0001));
    REQUIRE(__F_EQ_IN(oo_r(2, 1).imag(), 0, 0.0001));

    save_o = o;
    q2.inv(3, 3, o.data(), o2_r.data());
    REQUIRE(__F_EQ_IN(o2_r(0, 1), 1 * 9, 0.0001));
    REQUIRE(__F_EQ_IN(o2_r(1, 0), 3 * 9, 0.0001));

    q2.inv(3, 3, save_o.data(), o2_r.data());
    REQUIRE(__F_EQ_IN(o2_r(0, 1), 1 * 9, 0.0001));
    REQUIRE(__F_EQ_IN(o2_r(1, 0), 3 * 9, 0.0001));

    // r2r 2d
    fftw3::plan<double> p3, q3;
    Matrix<double, 3, 3, RowMajor> i_rr(3, 3), o_rr(3, 3), o2_rr(3, 3);

    i_rr << 0, 1, 2, 3, 4, 5, 6, 7, 8;
    p3.fwd<DCT_II, DCT_II>(3, 3, i_rr.data(), o_rr.data());
    // octave scaled values
    // cout << o_rr(0, 0) * std::sqrt(1.0/9)/4<< endl;
    // cout << o_rr(0, 1) * std::sqrt(2.0/9)/4<< endl;
    // cout << o_rr(1, 0) * std::sqrt(2.0/9)/4<< endl;
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 144, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -20.7846, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(1, 0), -62.3538, 0.0001));

    p3.fwd<DCT_II, DCT_II>(3, 3, i_rr.data(), o_rr.data());
    REQUIRE(__F_EQ_IN(o_rr(0, 0), 144, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(0, 1), -20.7846, 0.0001));
    REQUIRE(__F_EQ_IN(o_rr(1, 0), -62.3538, 0.0001));

    q3.inv<DCT_II, DCT_II>(3, 3, o_rr.data(), o2_rr.data());
    // scaled by 2n*2n, where n is 3
    REQUIRE(__F_EQ_IN(o2_rr(0, 1), 1 * 36, 0.0001));
    REQUIRE(__F_EQ_IN(o2_rr(1, 0), 3 * 36, 0.0001));
    REQUIRE(__F_EQ_IN(o2_rr(2, 2), 8 * 36, 0.0001));

    q3.inv<DCT_II, DCT_II>(3, 3, o_rr.data(), o2_rr.data());
    REQUIRE(__F_EQ_IN(o2_rr(0, 1), 1 * 36, 0.0001));
    REQUIRE(__F_EQ_IN(o2_rr(1, 0), 3 * 36, 0.0001));
    REQUIRE(__F_EQ_IN(o2_rr(2, 2), 8 * 36, 0.0001));
}
