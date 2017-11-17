#include <../test/test_util.h>
#include <Eigen/Dense>
#include <catch.hpp>
#include <iostream>
#include <polynomial/evaluation.h>

TEST_CASE("poly_eval")
{
    iexp::VectorXd v(3);
    v << 1.0, 0.5, 0.3;
    double x = 0.5;
    double ans = iexp::poly::eval(v.array(), x);
    REQUIRE(__D_EQ(1.0 + 0.5 * x + 0.3 * x * x, ans));

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

    // test compile
    iexp::Array<double, iexp::Dynamic, 1> a(10), b(10);
    b = iexp::poly::eval_deriv(a, x, 3);
    b = iexp::poly::eval_deriv(a + a + a, x, 3);
}
