#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <poly/evaluation.h>

TEST_CASE("poly_eval")
{
    iexp::VectorXd v(3);
    v << 1.0, 0.5, 0.3;
    double x = 0.5;
    double ans = iexp::poly::eval(v.array(), x);
    REQUIRE(__D_EQ(1.0 + 0.5 * x + 0.3 * x * x, ans));

    ans = iexp::poly::eval(v, x);
    REQUIRE(__D_EQ(1.0 + 0.5 * x + 0.3 * x * x, ans));

    iexp::Matrix<double, 2, 2> mm;
    mm << 1, 2, 3, 4;
    ans = iexp::poly::eval(mm, x);
    REQUIRE(__D_EQ(1 + 3 * x + 2 * x * x + 4 * x * x * x, ans));

    iexp::Vector4cd v2(4);
    v2[0] = std::complex<double>(-2.31, 0.44);
    v2[1] = std::complex<double>(4.21, -3.19);
    v2[2] = std::complex<double>(0.93, 1.04);
    v2[3] = std::complex<double>(-0.42, 0.68);
    std::complex<double> y(0.49, 0.95);
    std::complex<double> ans2 = iexp::poly::eval(v2.array(), y);
    REQUIRE(__D_EQ_IN(ans2.real(), 1.82462012, 0.00000001));
    REQUIRE(__D_EQ_IN(ans2.imag(), 2.30389412, 0.00000001));

    // test compile
    iexp::Array<double, 1, 10> a(10);
    ans = iexp::poly::eval(a.array(), x);
    ans = iexp::poly::eval((a.array() + a.array() + a.array()), x);

    // eval more
    iexp::Matrix<double, 3, 1> cc;
    cc << 1, 0, 3; // 1 + x^2
    iexp::Matrix<double, 2, 3> mx, mr;
    mx << 1, 2, 3, 4, 5, 6;
    mr = iexp::poly::eval(cc, mx);
    REQUIRE(__D_EQ(mr(0, 0), 4));
    REQUIRE(__D_EQ(mr(1, 0), 49));
    REQUIRE(__D_EQ(mr(1, 2), 109));

    iexp::Matrix<double, 2, 3, iexp::RowMajor> mx2 = mx, mr2;
    mr2 = iexp::poly::eval(cc, mx2);
    REQUIRE(__D_EQ(mr2(0, 0), 4));
    REQUIRE(__D_EQ(mr2(1, 0), 49));
    REQUIRE(__D_EQ(mr2(1, 2), 109));
}

TEST_CASE("poly_eval_derive")
{
    iexp::VectorXd v(6), v2(6);
    v << 1.0, -2.0, +3.0, -4.0, +5.0, -6.0;
    double x = -0.5;
    v2 = iexp::poly::eval_deriv(v.array(), x, 5);
    double y = v[0] + v[1] * x + v[2] * x * x + v[3] * x * x * x +
               v[4] * x * x * x * x + v[5] * x * x * x * x * x;
    REQUIRE(__D_EQ_IN(v2[0], y, 0.001));
    y = v[1] + 2.0 * v[2] * x + 3.0 * v[3] * x * x + 4.0 * v[4] * x * x * x +
        5.0 * v[5] * x * x * x * x;
    REQUIRE(__D_EQ_IN(v2[1], y, 0.001));
    y = 2.0 * v[2] + 3.0 * 2.0 * v[3] * x + 4.0 * 3.0 * v[4] * x * x +
        5.0 * 4.0 * v[5] * x * x * x;
    REQUIRE(__D_EQ_IN(v2[2], y, 0.001));
    y = 3.0 * 2.0 * v[3] + 4.0 * 3.0 * 2.0 * v[4] * x +
        5.0 * 4.0 * 3.0 * v[5] * x * x;
    REQUIRE(__D_EQ_IN(v2[3], y, 0.001));
    y = 4.0 * 3.0 * 2.0 * v[4] + 5.0 * 4.0 * 3.0 * 2.0 * v[5] * x;
    REQUIRE(__D_EQ_IN(v2[4], y, 0.001));
    y = 5.0 * 4.0 * 3.0 * 2.0 * v[5];
    REQUIRE(__D_EQ_IN(v2[5], y, 0.001));

    // matrix
    v2 = iexp::poly::eval_deriv(v, x, 5);
    REQUIRE(v2.size() == 6);
    y = v[0] + v[1] * x + v[2] * x * x + v[3] * x * x * x +
        v[4] * x * x * x * x + v[5] * x * x * x * x * x;
    REQUIRE(__D_EQ_IN(v2[0], y, 0.001));
    y = 5.0 * 4.0 * 3.0 * 2.0 * v[5];
    REQUIRE(__D_EQ_IN(v2[5], y, 0.001));

    // test compile
    iexp::Array<double, iexp::Dynamic, 1> a(10), b(10);
    b = iexp::poly::eval_deriv(a, x, 3);
    b = iexp::poly::eval_deriv(a + a + a, x, 3);

    // eval more
    iexp::Matrix<double, 3, 1> cc;
    cc << 1, 0, 3; // 1 + 3 * x^2
    iexp::Matrix<double, 2, 3> mx, mr;
    mx << 1, 2, 3, 4, 5, 6;
    mr = iexp::poly::eval(cc, mx);
    REQUIRE(__D_EQ(mr(0, 0), 4));
    REQUIRE(__D_EQ(mr(1, 0), 49));
    REQUIRE(__D_EQ(mr(1, 2), 109));

    // eval more derive
    iexp::Matrix<double, 3, 6> md1;
    md1 = iexp::poly::eval_deriv(cc, mx, 2);
    // eval 1 + 3*x^2, 6*x, 6 of [1,4,2,5,3,6]
    REQUIRE(__D_EQ(md1(0, 0), 4));
    REQUIRE(__D_EQ(md1(2, 0), 6));
    REQUIRE(__D_EQ(md1(0, 1), 49));
    REQUIRE(__D_EQ(md1(1, 1), 24));
    REQUIRE(__D_EQ(md1(0, 5), 109));
    REQUIRE(__D_EQ(md1(2, 5), 6));

    iexp::Matrix<double, 6, 3> md2;
    md2 = iexp::poly::eval_deriv<true>(cc, mx, 2);
    REQUIRE(__D_EQ(md2(0, 0), 4));
    REQUIRE(__D_EQ(md2(0, 2), 6));
    REQUIRE(__D_EQ(md2(1, 0), 49));
    REQUIRE(__D_EQ(md2(1, 1), 24));
    REQUIRE(__D_EQ(md2(5, 0), 109));
    REQUIRE(__D_EQ(md2(5, 2), 6));
}
