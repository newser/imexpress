#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <poly/divided_difference.h>
#include <poly/evaluation.h>

TEST_CASE("poly_dd")
{
    iexp::VectorXd xa(7), ya(7), dd(7);
    xa << 0.16, 0.97, 1.94, 2.74, 3.58, 3.73, 4.70;
    ya << 0.73, 1.11, 1.49, 1.84, 2.30, 2.41, 3.07;
    dd = iexp::poly::dd(xa.array(), ya.array());
    REQUIRE(dd.size() == 7);
    REQUIRE(__D_EQ_IN(dd[0], 7.30000000000000e-01, 0.000000000000001));
    REQUIRE(__D_EQ_IN(dd[1], 4.69135802469136e-01, 0.000000000000001));
    REQUIRE(__D_EQ_IN(dd[2], -4.34737219941284e-02, 0.000000000000001));
    REQUIRE(__D_EQ_IN(dd[3], 2.68681098870099e-02, 0.000000000000001));
    REQUIRE(__D_EQ_IN(dd[4], -3.22937056934996e-03, 0.000000000000001));
    REQUIRE(__D_EQ_IN(dd[5], 6.12763259971375e-03, 0.000000000000001));
    REQUIRE(__D_EQ_IN(dd[6], -6.45402453527083e-03, 0.000000000000001));

    dd = iexp::poly::dd(xa, ya);
    REQUIRE(dd.size() == 7);
    REQUIRE(__D_EQ_IN(dd[0], 7.30000000000000e-01, 0.000000000000001));
    REQUIRE(__D_EQ_IN(dd[6], -6.45402453527083e-03, 0.000000000000001));

    iexp::Matrix<double, 2, 2, iexp::RowMajor> xxa, yya, ddd;
    xxa << 0.16, 0.97, 1.94, 2.74;
    yya << 0.73, 1.11, 1.49, 1.84;
    ddd = iexp::poly::dd(xxa, yya);
    REQUIRE(__D_EQ_IN(ddd(0, 0), 7.30000000000000e-01, 0.000000000000001));
    REQUIRE(__D_EQ_IN(ddd(0, 1), 4.69135802469136e-01, 0.000000000000001));
    REQUIRE(__D_EQ_IN(ddd(1, 0), -4.34737219941284e-02, 0.000000000000001));
    REQUIRE(__D_EQ_IN(ddd(1, 1), 2.68681098870099e-02, 0.000000000000001));

    // dd eval
    for (int i = 0; i < xa.size(); ++i) {
        double ans = iexp::poly::dd_eval(dd.array(), xa.array(), xa[i]);
        REQUIRE(__D_EQ_IN(ya[i], ans, 0.000000000000001));
    }

    iexp::MatrixXd x2(2, 3), y2(2, 3), ans2(2, 3);
    x2 << 0.16, 0.97, 1.94, 2.74, 3.58, 3.73;
    y2 << 0.73, 1.11, 1.49, 1.84, 2.30, 2.41;
    ans2 = iexp::poly::dd_eval(dd, xa, x2);
    REQUIRE(__D_EQ9(ans2(0, 0), y2(0, 0)));
    REQUIRE(__D_EQ9(ans2(0, 2), y2(0, 2)));
    REQUIRE(__D_EQ9(ans2(1, 0), y2(1, 0)));
    REQUIRE(__D_EQ9(ans2(1, 2), y2(1, 2)));

    // dd taylor
    iexp::VectorXd c = iexp::poly::dd_taylor(1.5, dd.array(), xa.array());
    for (int i = 0; i < xa.size(); ++i) {
        double ans = iexp::poly::eval(c.array(), xa[i] - 1.5);
        REQUIRE(__D_EQ_IN(ans, ya[i], 1e-10));
    }

    // compile
    dd = iexp::poly::dd((xa.array() + ya.array()) * ya.array(),
                        (xa.array() + ya.array()) * ya.array());
    double ans = iexp::poly::dd_eval((xa.array() + ya.array()) * ya.array(),
                                     (xa.array() + ya.array()) * ya.array(),
                                     0);
}

TEST_CASE("poly_dd_hermit")
{
    iexp::ArrayXd xa(3), ya(3), dya(3), ans(6);
    xa << 1.3, 1.6, 1.9;
    ya << 0.6200860, 0.4554022, 0.2818186;
    dya << -0.5220232, -0.5698959, -0.5811571;
    ans << 6.200860000000e-01, -5.220232000000e-01, -8.974266666667e-02,
        6.636555555556e-02, 2.666666666662e-03, -2.774691357989e-03;

    iexp::ArrayXd h = iexp::poly::dd_hermit(xa, ya, dya);
    REQUIRE(h.size() == 6);
    for (int i = 0; i < h.size(); ++i) {
        REQUIRE(__D_EQ_IN(h[i], ans[i], 1e-10));
    }

    iexp::RowVectorXd xa2, ya2, dya2;
    xa2 = xa, ya2 = ya, dya2 = dya;
    iexp::RowVectorXd h2 = iexp::poly::dd_hermit(xa2, ya2, dya2);
    REQUIRE(h2.size() == 6);
    for (int i = 0; i < h2.size(); ++i) {
        REQUIRE(__D_EQ_IN(h2[i], ans[i], 1e-10));
    }

    // compile
    h = iexp::poly::dd_hermit(xa * ya + dya, xa * ya + dya, xa * ya + dya);
}
