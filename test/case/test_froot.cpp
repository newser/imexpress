#include <catch.hpp>
#include <common/gslvec.h>
#include <froot/dfroot.h>
#include <froot/froot.h>
#include <froot/mdfroot.h>
#include <froot/mfroot.h>
#include <iostream>
#include <math/constant.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_froot")
{
    // 0 root
    for (int i = 0; i < 3; ++i) {
        try {
            // root: 1, 2
            froot
            fr([](double x)
                   -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
               0.1,
               0.99,
               (froot::type)i);
            REQUIRE(false);
        } catch (...) {
        }
    }

    // 1 root
    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        froot fr([](double x)
                     -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                 0.1,
                 1.2,
                 (froot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE(__D_EQ6(rt, 1));

        rt = fr.find<false>(0, 1e-6);
        REQUIRE(__D_EQ6(rt, 1));

        rt = fr.find(1e-6);
        REQUIRE(__D_EQ6(rt, 1));
    }

    // 2 root
    for (int i = 1; i < 3; ++i) {
        try {
            // root: 1, 2
            froot
            fr([](double x)
                   -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
               0.1,
               2.2,
               (froot::type)i);
            REQUIRE(false);
        } catch (...) {
        }
    }

    // multiple root in interval
    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        froot fr([](double x)
                     -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                 0.1,
                 3.2,
                 (froot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));

        rt = fr.find<false>(0, 1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));

        rt = fr.find(1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));
    }
}

TEST_CASE("test_dfroot")
{
    // 1 root
    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        dfroot fr([](double x)
                      -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                  [](double x) -> double { return 3 * x * x - 12 * x + 11; },
                  [](double x, double &f, double &df) -> void {
                      f = x * x * x - 6 * x * x + 11 * x - 6;
                      df = 3 * x * x - 12 * x + 11;
                  },
                  -0.1,
                  (dfroot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE(__D_EQ6(rt, 1));

        rt = fr.find(1e-6);
        REQUIRE(__D_EQ6(rt, 1));
    }

    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        dfroot fr([](double x)
                      -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                  [](double x) -> double { return 3 * x * x - 12 * x + 11; },
                  [](double x, double &f, double &df) -> void {
                      f = x * x * x - 6 * x * x + 11 * x - 6;
                      df = 3 * x * x - 12 * x + 11;
                  },
                  3.2,
                  (dfroot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));

        rt = fr.find(1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));
    }
}

TEST_CASE("test_gslvec")
{
    {
        VectorXcd v(3);

        // noncopy
        gslvec<std::complex<double>> gv(v, false);
        gsl_vector_complex *pgv = gv.gsl_vector();
        v[0] = std::complex<double>(1, 2);
        gsl_complex vv = gsl_vector_complex_get(pgv, 0);
        REQUIRE(vv.dat[0] == 1);
        REQUIRE(vv.dat[1] == 2);
        vv.dat[0] = 3;
        vv.dat[1] = 4;
        gsl_vector_complex_set(pgv, 1, vv);
        REQUIRE(v[1] == std::complex<double>(3, 4));

        // copy
        v.setConstant(std::complex<double>(12, 34));
        gslvec<std::complex<double>> gv2(v, true);
        pgv = gv2.gsl_vector();
        v[0] = std::complex<double>(23, 45);
        vv = gsl_vector_complex_get(pgv, 0);
        REQUIRE(vv.dat[0] == 12);
        REQUIRE(vv.dat[1] == 34);
        vv.dat[0] = 34;
        vv.dat[1] = 45;
        gsl_vector_complex_set(pgv, 1, vv);
        REQUIRE(v[1] == std::complex<double>(12, 34));

        // can not copy
        v.setConstant(std::complex<double>(12, 34));
        gslvec<std::complex<double>> gv3(v + v);
        pgv = gv3.gsl_vector();
        v[0] = std::complex<double>(23, 45);
        vv = gsl_vector_complex_get(pgv, 0);
        REQUIRE(vv.dat[0] == 24);
        REQUIRE(vv.dat[1] == 68);
    }

    {
        VectorXd v(3);

        gslvec<double> gv(v, false);
        gsl_vector *pgv = gv.gsl_vector();
        v[0] = 1;
        REQUIRE(gsl_vector_get(pgv, 0) == 1);
        gsl_vector_set(pgv, 1, 9);
        REQUIRE(v[1] == 9);
    }
}

TEST_CASE("test_mfroot")
{
    for (int i = 0; i < 4; ++i) {
        Vector2d x;
        x << -10.0, -5.0;

        mfroot fr(
            [](Map<const VectorXd> &x, Map<VectorXd> &f) -> bool {
                f(0) = 1 - x(0);
                f(1) = 10 * (x(1) - x(0) * x(0));
                return true;
            },
            x,
            (mfroot::type)i);

        Vector2d rt = fr.find(0, 1e-6);
        REQUIRE(__D_EQ6(rt(0), 1));
        REQUIRE(__D_EQ6(rt(1), 1));

        mfroot fr2(
            [](Map<const VectorXd> &x, Map<VectorXd> &f) -> bool {
                f(0) = 1 - x(0);
                f(1) = 10 * (x(1) - x(0) * x(0));
                return true;
            },
            x,
            (mfroot::type)i);

        rt = fr2.find(1e-6);
        REQUIRE(__D_EQ6(rt(0), 1));
        REQUIRE(__D_EQ6(rt(1), 1));
    }
}

bool jac(Map<const VectorXd> &x, Map<RowMatrixXd> &jac)
{
    jac(0, 0) = -1;
    jac(0, 1) = 0;
    jac(1, 0) = -20 * x(0);
    jac(1, 1) = 10;
    return true;
}

bool fdf(Map<const VectorXd> &x, Map<VectorXd> &f, Map<RowMatrixXd> &jac)
{
    f(0) = 1 - x(0);
    f(1) = 10 * (x(1) - x(0) * x(0));

    jac(0, 0) = -1;
    jac(0, 1) = 0;
    jac(1, 0) = -20 * x(0);
    jac(1, 1) = 10;
    return true;
}

TEST_CASE("test_mdfroot")
{
    for (int i = 0; i < 4; ++i) {
        Vector2d x;
        x << -10.0, -5.0;

        mdfroot fr(
            [](Map<const VectorXd> &x, Map<VectorXd> &f) -> bool {
                f(0) = 1 - x(0);
                f(1) = 10 * (x(1) - x(0) * x(0));
                return true;
            },
            jac,
            fdf,
            x,
            (mdfroot::type)i);

        Vector2d rt = fr.find(0, 1e-6);
        REQUIRE(__D_EQ6(rt(0), 1));
        REQUIRE(__D_EQ6(rt(1), 1));

        mfroot fr2(
            [](Map<const VectorXd> &x, Map<VectorXd> &f) -> bool {
                f(0) = 1 - x(0);
                f(1) = 10 * (x(1) - x(0) * x(0));
                return true;
            },
            x,
            (mfroot::type)i);

        rt = fr2.find(1e-6);
        REQUIRE(__D_EQ6(rt(0), 1));
        REQUIRE(__D_EQ6(rt(1), 1));
    }
}
