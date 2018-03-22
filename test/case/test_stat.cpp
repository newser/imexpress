#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <stats/mean.h>
#include <stats/std.h>
#include <stats/var.h>

TEST_CASE("stat_mean")
{
    double v, g;

    iexp::ArrayXd c;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::mean(c);
    g = gsl_stats_mean(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::mean(c2);
    g = gsl_stats_int_mean(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::mean(c + c2.cast<double>());
}

TEST_CASE("stat_var")
{
    double v, g;

    iexp::ArrayXd c;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::var(c);
    g = gsl_stats_variance(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::var(c, 1);
    g = gsl_stats_variance_m(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::var(c2);
    g = gsl_stats_int_variance(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::var(c2, 2);
    g = gsl_stats_int_variance_m(c2.data(), 1, c2.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::var(c + c2.cast<double>());
    v = iexp::stats::var(c + c2.cast<double>(), 2);
}

TEST_CASE("stat_std")
{
    double v, g;

    iexp::ArrayXd c;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::std(c);
    g = gsl_stats_sd(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::std(c, 1);
    g = gsl_stats_sd_m(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::std(c2);
    g = gsl_stats_int_sd(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::std(c2, 2);
    g = gsl_stats_int_sd_m(c2.data(), 1, c2.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::std(c + c2.cast<double>());
}
