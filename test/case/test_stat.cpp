#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <stats/mean.h>

TEST_CASE("stat_common")
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
