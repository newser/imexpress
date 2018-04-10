#include <catch.hpp>
#include <integral/cquad.h>
#include <integral/miser.h>
#include <integral/monte.h>
#include <integral/qag.h>
#include <integral/qagi.h>
#include <integral/qagp.h>
#include <integral/qags.h>
#include <integral/qawc.h>
#include <integral/qawf.h>
#include <integral/qawo.h>
#include <integral/qaws.h>
#include <integral/qng.h>
#include <integral/vegas.h>
#include <iostream>
#include <test_util.h>

using namespace iexp;
using namespace std;

static double foo(double x)
{
    return x * x;
}

TEST_CASE("integral_func")
{
    integral::unary_func<double> f([](double x) -> double { return 2 * x; });
    int dummy[10]{1, 2, 3, 4, 5, 6, 7, 8, 9};
    double v = f.gsl()->function(1.2, f.gsl()->params);
    REQUIRE(v == 2.4);

    // integral::unary_func<double> f2([](float x) -> float { return x; });

    class t_op
    {
      public:
        double operator()(double x)
        {
            return x + 2;
        }
    };
    integral::unary_func<double> f2((t_op()));
    int dummy2[10]{1, 2, 3, 4, 5, 6, 7, 8, 9};
    v = f2.gsl()->function(1.2, f2.gsl()->params);
    REQUIRE(v == 3.2);

    integral::unary_func<double> f3(foo);
    int dummy3[10]{1, 2, 3, 4, 5, 6, 7, 8, 9};
    v = f3.gsl()->function(1.2, f3.gsl()->params);
    REQUIRE(v == 1.44);
}

TEST_CASE("integral_qng")
{
    double result, abserr;
    size_t neval;

    integral::qng q(1e-1, 0.0);
    result = q([](double x) { return pow(x, 2.6) * log(1 / x); },
               0,
               1,
               &abserr,
               &neval);
    REQUIRE(__D_EQ9(result, 7.716049379303083211E-02));
    REQUIRE(__D_EQ9(abserr, 9.424302199601294244E-08));
    REQUIRE(neval == 21);

    REQUIRE(q.epsabs() == 1e-1);
    REQUIRE(q.epsrel() == 0.0);

    result = q([](double x) { return pow(x, 2.6) * log(1 / x); },
               1,
               0,
               &abserr,
               &neval);
    REQUIRE(__D_EQ9(result, -7.716049379303083211E-02));
    REQUIRE(__D_EQ9(abserr, 9.424302199601294244E-08));
    REQUIRE(neval == 21);

    result = q([](double x) { return pow(x, 2.6) * log(1 / x); }, 0, 1);
    REQUIRE(__D_EQ9(result, 7.716049379303083211E-02));
}

TEST_CASE("integral_qag")
{
    double result, abserr;

    integral::qag q(0.0, 1e-10, integral::key::GAUSS15, 1000);
    result =
        q([](double x) { return pow(x, 2.6) * log(1 / x); }, 0.0, 1.0, &abserr);
    REQUIRE(__D_EQ_IN(result, 7.716049382715854665E-02, 1E-15));
    REQUIRE(__D_EQ_IN(abserr, 6.679384885865053037E-12, 1E-6));

    REQUIRE(q.epsabs() == 0.0);
    REQUIRE(q.epsrel() == 1e-10);
    REQUIRE(q.key() == integral::key::GAUSS15);

    result =
        q([](double x) { return pow(x, 2.6) * log(1 / x); }, 1.0, 0.0, &abserr);
    REQUIRE(__D_EQ_IN(result, -7.716049382715854665E-02, 1E-15));
    REQUIRE(__D_EQ_IN(abserr, 6.679384885865053037E-12, 1E-6));

    result = q([](double x) { return pow(x, 2.6) * log(1 / x); }, 1.0, 0.0);
    REQUIRE(__D_EQ_IN(result, -7.716049382715854665E-02, 1E-15));
}

TEST_CASE("integral_qags")
{
    double result, abserr;

    integral::qags q(0.0, 1e-10, 1000);
    result =
        q([](double x) { return pow(x, 2.6) * log(1 / x); }, 0.0, 1.0, &abserr);
    REQUIRE(__D_EQ_IN(result, 7.716049382715854665E-02, 1E-15));
    REQUIRE(__D_EQ_IN(abserr, 6.679384885865053037E-12, 1E-6));

    REQUIRE(q.epsabs() == 0.0);
    REQUIRE(q.epsrel() == 1e-10);

    result =
        q([](double x) { return pow(x, 2.6) * log(1 / x); }, 1.0, 0.0, &abserr);
    REQUIRE(__D_EQ_IN(result, -7.716049382715854665E-02, 1E-15));
    REQUIRE(__D_EQ_IN(abserr, 6.679384885865053037E-12, 1E-6));
}

TEST_CASE("integral_qagp")
{
    double pts[4] = {0.0, 1.0, sqrt(2.0), 3.0};
    double result, abserr;

    integral::qagp q(0.0, 1.0e-3, 1000);
    result = q(
        [](double x) {
            double x2 = x * x;
            double x3 = x * x2;
            return x3 * log(fabs((x2 - 1.0) * (x2 - 2.0)));
        },
        pts,
        4,
        &abserr);
    REQUIRE(__D_EQ_IN(result, 5.274080611672716401E+01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.755703848687062418E-04, 1E-5));

    REQUIRE(q.epsabs() == 0.0);
    REQUIRE(q.epsrel() == 1.0e-3);

    result = q(
        [](double x) {
            double x2 = x * x;
            double x3 = x * x2;
            return x3 * log(fabs((x2 - 1.0) * (x2 - 2.0)));
        },
        {0.0, 1.0, sqrt(2.0), 3.0},
        &abserr);
    REQUIRE(__D_EQ_IN(result, 5.274080611672716401E+01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.755703848687062418E-04, 1E-5));
}

TEST_CASE("integral_qagi")
{
    double result, abserr;

    integral::qagi q(1.0e-7, 0.0, 1000);
    result = q([](double x) { return exp(-x - x * x); }, &abserr);
    REQUIRE(__D_EQ_IN(result, 2.275875794468747770E+00, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 7.436490118267390744E-09, 1E-5));

    REQUIRE(q.epsabs() == 1.0e-7);
    REQUIRE(q.epsrel() == 0.0);

    result = q([](double x) { return exp(-x - x * x); }, &abserr);
    REQUIRE(__D_EQ_IN(result, 2.275875794468747770E+00, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 7.436490118267390744E-09, 1E-5));

    result = q([](double x) { return exp(-x - x * x); });
    REQUIRE(__D_EQ_IN(result, 2.275875794468747770E+00, 1E-14));
}

TEST_CASE("integral_qagil")
{
    double result, abserr;

    integral::qagil q(1.0e-7, 0.0, 1000);
    result = q([](double x) { return exp(x); }, 1.0, &abserr);
    REQUIRE(__D_EQ_IN(result, 2.718281828459044647E+00, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.588185109253204805E-10, 1E-5));

    REQUIRE(q.epsabs() == 1.0e-7);
    REQUIRE(q.epsrel() == 0.0);

    result = q([](double x) { return exp(x); }, 1.0, &abserr);
    REQUIRE(__D_EQ_IN(result, 2.718281828459044647E+00, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.588185109253204805E-10, 1E-5));
}

TEST_CASE("integral_qagiu")
{
    double result, abserr;

    integral::qagiu q(0.0, 1.0e-3, 1000);
    result = q([](double x) { return log(x) / (1.0 + 100.0 * x * x); },
               0.0,
               &abserr);
    REQUIRE(__D_EQ_IN(result, -3.616892186127022568E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 3.016716913328831851E-06, 1E-5));

    REQUIRE(q.epsabs() == 0.0);
    REQUIRE(q.epsrel() == 1.0e-3);

    result = q([](double x) { return log(x) / (1.0 + 100.0 * x * x); },
               0.0,
               &abserr);
    REQUIRE(__D_EQ_IN(result, -3.616892186127022568E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 3.016716913328831851E-06, 1E-5));
}

TEST_CASE("integral_qawc")
{
    double result, abserr;

    integral::qawc q(0.0, 1.0e-3, 1000);
    result = q([](double x) { return 1.0 / (5.0 * x * x * x + 6.0); },
               -1.0,
               5.0,
               0.0,
               &abserr);
    REQUIRE(__D_EQ_IN(result, -8.994400695837000137E-02, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.185290176227023727E-06, 1E-6));

    REQUIRE(q.epsabs() == 0.0);
    REQUIRE(q.epsrel() == 1.0e-3);

    result = q([](double x) { return 1.0 / (5.0 * x * x * x + 6.0); },
               -1.0,
               5.0,
               0.0,
               &abserr);
    REQUIRE(__D_EQ_IN(result, -8.994400695837000137E-02, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.185290176227023727E-06, 1E-6));
}

TEST_CASE("integral_qaws")
{
    double result, abserr;

    integral::qaws_table t(0.0, 0.0, 1, 0);

    integral::qaws q(

        0.0, 1.0e-7, 1000);
    result = q(
        [](double x) -> double {
            if (x == 0.0) {
                return 0.0;
            } else {
                double u = log(x);
                double v = 1 + u * u;

                return 1.0 / (v * v);
            }
        },
        t,
        0.0,
        1.0,
        &abserr);
    REQUIRE(__D_EQ_IN(result, -1.892751853489401670E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.129133712015747658E-08, 1E-6));

    t.alpha(-0.5).beta(-0.3).mu(1).nu(1);
    t.reload();

    q.epsabs(0.0);
    q.epsrel(1.0e-7);
    result = q(
        [](double x) -> double {
            if (x == 0.0) {
                return 0.0;
            } else {
                double u = log(x);
                double v = 1 + u * u;

                return 1.0 / (v * v);
            }
        },
        t,
        0.0,
        1.0,
        &abserr);
    REQUIRE(__D_EQ_IN(result, 3.159922862811048172E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 2.336183482198144595E-08, 1E-6));
}

TEST_CASE("integral_qawo_table")
{
    integral::qawo_table t(10.0 * M_PI, 1.0, true, 1000);
    t.reload();
    REQUIRE(t.gsl() != nullptr);

    t.omega(9.0 * M_PI).length(2.0).sine(false);
    t.reload();
    REQUIRE(t.gsl() != nullptr);
}

TEST_CASE("integral_qawo")
{
    double result, abserr;

    integral::qawo_table t(10.0 * M_PI, 1.0, true, 1000);

    integral::qawo q(0.0, 1e-7, 1000);
    result = q(
        [](double x) -> double {
            if (x == 0.0) {
                return 0.0;
            }
            return log(x);
        },
        t,
        0.0,
        &abserr);
    REQUIRE(__D_EQ_IN(result, -1.281368483991674190E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 6.875028324415666248E-12, 1E-3));

    t.length(-1.0);
    t.reload();

    q.epsabs(0.0).epsrel(1.0e-7);
    result = q(
        [](double x) -> double {
            if (x == 0.0) {
                return 0.0;
            }
            return log(x);
        },
        t,
        1.0,
        &abserr);
    REQUIRE(__D_EQ_IN(result, 1.281368483991674190E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 6.875028324415666248E-12, 1E-3));
}

TEST_CASE("integral_qawf")
{
    double result, abserr;

    integral::qawo_table t(M_PI / 2.0, 1.0, false, 1000);

    integral::qawf q(

        1e-7, 1000);
    result = q(
        [](double x) -> double {
            if (x == 0.0) {
                return 0;
            }
            return 1 / sqrt(x);
        },
        t,
        0.0,
        &abserr);
    REQUIRE(__D_EQ_IN(result, 9.999999999279802765E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.556289974669056164E-08, 1E-3));

    result = q(
        [](double x) -> double {
            if (x == 0.0) {
                return 0;
            }
            return 1 / sqrt(x);
        },
        t,
        0.0,
        &abserr);
    REQUIRE(__D_EQ_IN(result, 9.999999999279802765E-01, 1E-14));
    REQUIRE(__D_EQ_IN(abserr, 1.556289974669056164E-08, 1E-3));
}

static double foo_1d(double x)
{
    return x * 3;
}

static double foo_2d(double *x, size_t d)
{
    REQUIRE(d == 2);
    return x[0] + 3 * x[1];
}

TEST_CASE("integral_monte_func")
{
    ////////////////////////////
    // 1d
    ////////////////////////////

    integral::monte_func<double> f_1d([](double x) -> double { return 5 * x; });
    double d[2]{9, 8};
    double v = f_1d.gsl()->f(d, 1, f_1d.gsl()->params);
    REQUIRE(v == 45);

    class t_op2
    {
      public:
        double operator()(double x)
        {
            return 4 * x;
        }
    };
    integral::monte_func<double> f2_1d((t_op2()));
    v = f2_1d.gsl()->f(d, 2, f2_1d.gsl()->params);
    REQUIRE(v == 36);

    integral::monte_func<double> f3_1d(foo_1d);
    v = f3_1d.gsl()->f(d, 2, f3_1d.gsl()->params);
    REQUIRE(v == 27);

    ////////////////////////////
    // 2d
    ////////////////////////////

    integral::monte_func<double> f(2, [](double *x, size_t d) -> double {
        REQUIRE(d == 2);
        return 3 * x[0] * x[0] + 2 * x[0] * x[1] + x[1] * x[1];
    });
    v = f.gsl()->f(d, 2, f.gsl()->params);
    REQUIRE(v == 451);

    // integral::unary_func<double> f2([](float x) -> float { return x; });

    class t_op
    {
      public:
        double operator()(double *x, size_t d)
        {
            REQUIRE(d == 2);
            return x[0] + 2 * x[1];
        }
    };
    integral::monte_func<double> f2(2, t_op());
    v = f2.gsl()->f(d, 2, f2.gsl()->params);
    REQUIRE(v == 25);

    integral::monte_func<double> f3(2, foo_2d);
    v = f3.gsl()->f(d, 2, f3.gsl()->params);
    REQUIRE(v == 33);
}

TEST_CASE("integral_monte")
{
    double xl[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double xu[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    // 2d

    integral::monte m(2);
    double e;
    double v = m(
        [](double *x, size_t d) -> double {
            double prod = 1.0;
            unsigned int i;

            for (i = 0; i < d; ++i) {
                prod *= 2.0 * x[i];
            }

            return prod;
        },
        xl,
        xu,
        2,
        3333,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    m.reset();
    v = m([](double *x,
             size_t d) -> double { return (x[0] > 0.1 && x[0] < 0.9) ? 1 : 0; },
          xl,
          xu,
          2,
          100000,
          &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(__D_EQ_IN(e, 1.264e-3, 1e-3));

    // nd

    integral::monte m_n(4);
    v = m_n(
        [](double *x, size_t d) -> double {
            double prod = 1.0;
            unsigned int i;

            for (i = 0; i < d; ++i) {
                prod *= 2.0 * x[i];
            }

            return prod;
        },
        xl,
        xu,
        4,
        21604,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    v = m_n([](double *x, size_t d)
                -> double { return (x[0] > 0.1 && x[0] < 0.9) ? 1 : 0; },
            xl,
            xu,
            4,
            100000,
            &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(__D_EQ_IN(e, 1.264e-3, 1e-3));

    // 1d

    integral::monte m_1(1);
    v = m_1(
        [](double x) -> double {
            double prod = 1.0;

            prod *= 2.0 * x;

            return prod;
        },
        xl[0],
        xu[0],
        21604,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    v = m_1([](double x) -> double { return (x > 0.1 && x < 0.9) ? 1 : 0; },
            xl[0],
            xu[0],
            100000,
            &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(__D_EQ_IN(e, 1.264e-3, 1e-3));
}

TEST_CASE("integral_miser")
{
    double xl[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double xu[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    // 2d

    integral::miser m(2);
    double e;
    double v = m(
        [](double *x, size_t d) -> double {
            double prod = 1.0;
            unsigned int i;

            for (i = 0; i < d; ++i) {
                prod *= 2.0 * x[i];
            }

            return prod;
        },
        xl,
        xu,
        2,
        3333,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    m.reset();
    v = m([](double *x,
             size_t d) -> double { return (x[0] > 0.1 && x[0] < 0.9) ? 1 : 0; },
          xl,
          xu,
          2,
          100000,
          &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(e < 5 * 1.264e-3);

    // nd

    integral::miser m_n(4);
    v = m_n(
        [](double *x, size_t d) -> double {
            double prod = 1.0;
            unsigned int i;

            for (i = 0; i < d; ++i) {
                prod *= 2.0 * x[i];
            }

            return prod;
        },
        xl,
        xu,
        4,
        21604,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    v = m_n([](double *x, size_t d)
                -> double { return (x[0] > 0.1 && x[0] < 0.9) ? 1 : 0; },
            xl,
            xu,
            4,
            100000,
            &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(__D_EQ_IN(e, 1.264e-3, 1e-3));

    // 1d

    integral::miser m_1(1);
    v = m_1(
        [](double x) -> double {
            double prod = 1.0;

            prod *= 2.0 * x;

            return prod;
        },
        xl[0],
        xu[0],
        21604,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    v = m_1([](double x) -> double { return (x > 0.1 && x < 0.9) ? 1 : 0; },
            xl[0],
            xu[0],
            100000,
            &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(e < 5 * 1.264e-3);

    m_1.alpha(0.1).min_calls(10);
}

TEST_CASE("integral_vegas")
{
    double xl[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double xu[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    // 2d

    integral::vegas m(2);
    double e;
    double v = m(
        [](double *x, size_t d) -> double {
            double prod = 1.0;
            unsigned int i;

            for (i = 0; i < d; ++i) {
                prod *= 2.0 * x[i];
            }

            return prod;
        },
        xl,
        xu,
        2,
        3333,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    v = m([](double *x,
             size_t d) -> double { return (x[0] > 0.1 && x[0] < 0.9) ? 1 : 0; },
          xl,
          xu,
          2,
          100000,
          &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(e < 5 * 1.264e-3);

    // nd

    integral::vegas m_n(4);
    v = m_n(
        [](double *x, size_t d) -> double {
            double prod = 1.0;
            unsigned int i;

            for (i = 0; i < d; ++i) {
                prod *= 2.0 * x[i];
            }

            return prod;
        },
        xl,
        xu,
        4,
        21604,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    m_n.reset();
    v = m_n([](double *x, size_t d)
                -> double { return (x[0] > 0.1 && x[0] < 0.9) ? 1 : 0; },
            xl,
            xu,
            4,
            100000,
            &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(__D_EQ_IN(e, 1.264e-3, 1e-3));

    // 1d

    integral::vegas m_1(1);
    v = m_1(
        [](double x) -> double {
            double prod = 1.0;

            prod *= 2.0 * x;

            return prod;
        },
        xl[0],
        xu[0],
        21604,
        &e);
    REQUIRE(__D_EQ2(v, 1.0));
    REQUIRE(__D_EQ2(e, 0.01));

    v = m_1([](double x) -> double { return (x > 0.1 && x < 0.9) ? 1 : 0; },
            xl[0],
            xu[0],
            100000,
            &e);
    REQUIRE(__D_EQ_IN(v, 0.8, 0.1));
    REQUIRE(e < 5 * 1.264e-3);

    m_1.alpha(0.1).iterations(10);
}

TEST_CASE("integral_cquad")
{
    double result, abserr;
    size_t neval;

    integral::cquad q(0.0, 1e-12, 200);
    result = q([](double x) { return exp(x); }, 0.0, 1.0, &abserr, &neval);
    REQUIRE(__D_EQ_IN(result, 1.7182818284590452354, 1E-12));
    REQUIRE(fabs(result - 1.7182818284590452354) <= abserr);

    result = q([](double x) { return x >= 0.3; }, 0.0, 1.0, &abserr, &neval);
    REQUIRE(__D_EQ_IN(result, 0.7, 1E-12));
    REQUIRE(fabs(result - 0.7) <= abserr);
}
