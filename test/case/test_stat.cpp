#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <stats/autocorr.h>
#include <stats/corrcoef.h>
#include <stats/cov.h>
#include <stats/kurtosis.h>
#include <stats/mean.h>
#include <stats/skewness.h>
#include <stats/std.h>
#include <stats/tss.h>
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

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::uvar(c, 1);
    g = gsl_stats_variance_with_fixed_mean(c.data(), 1, c.size(), 1);
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

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::uvar(c2, 2);
    g = gsl_stats_int_variance_with_fixed_mean(c2.data(), 1, c2.size(), 2);
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

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::ustd(c, 1);
    g = gsl_stats_sd_with_fixed_mean(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::abstd(c);
    g = gsl_stats_absdev(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::abstd(c, 2);
    g = gsl_stats_absdev_m(c.data(), 1, c.size(), 2);
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

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::ustd(c2, 2);
    g = gsl_stats_int_sd_with_fixed_mean(c2.data(), 1, c2.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::abstd(c2);
    g = gsl_stats_int_absdev(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::abstd(c2, 2);
    g = gsl_stats_int_absdev_m(c2.data(), 1, c2.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::std(c + c2.cast<double>());
}

TEST_CASE("stat_tss")
{
    double v, g;

    iexp::ArrayXd c;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::tss(c);
    g = gsl_stats_tss(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::tss(c, 1);
    g = gsl_stats_tss_m(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::tss(c2);
    g = gsl_stats_int_sd(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::tss(c2, 2);
    g = gsl_stats_int_sd_m(c2.data(), 1, c2.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::tss(c + c2.cast<double>());
}

TEST_CASE("stat_skew")
{
    double v, g;

    iexp::ArrayXd c;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::skewness(c);
    g = gsl_stats_skew(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::skewness(c, 1, 2);
    g = gsl_stats_skew_m_sd(c.data(), 1, c.size(), 1, 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::skewness(c2);
    g = gsl_stats_int_skew(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::skewness(c2, 2, 3);
    g = gsl_stats_int_skew_m_sd(c2.data(), 1, c2.size(), 2, 3);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::skewness(c + c2.cast<double>());
}

TEST_CASE("stat_kurtosis")
{
    double v, g;

    iexp::ArrayXd c;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::kurtosis(c);
    g = gsl_stats_kurtosis(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::kurtosis(c, 1, 2);
    g = gsl_stats_kurtosis_m_sd(c.data(), 1, c.size(), 1, 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::kurtosis(c2);
    g = gsl_stats_int_kurtosis(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::kurtosis(c2, 2, 3);
    g = gsl_stats_int_kurtosis_m_sd(c2.data(), 1, c2.size(), 2, 3);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::kurtosis(c + c2.cast<double>());
}

TEST_CASE("stat_autocorr")
{
    double v, g;

    iexp::ArrayXd c;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::autocorr(c);
    g = gsl_stats_lag1_autocorrelation(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::autocorr(c, 1);
    g = gsl_stats_lag1_autocorrelation_m(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::autocorr(c2);
    g = gsl_stats_int_lag1_autocorrelation(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::autocorr(c2, 1);
    g = gsl_stats_int_lag1_autocorrelation_m(c2.data(), 1, c2.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::mean(c + c2.cast<double>());
}

TEST_CASE("stat_cov")
{
    double v, g;

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::cov(c, d);
    g = gsl_stats_covariance(c.data(), 1, d.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2, d2;
    c2 = iexp::ArrayXi::Random(10);
    d2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::cov(c2, d2);
    g = gsl_stats_int_covariance(c2.data(), 1, d2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::cov(c, d, 1, 2);
    g = gsl_stats_covariance_m(c.data(), 1, d.data(), 1, c.size(), 1, 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    d2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::cov(c2, d2, 3, 4);
    g = gsl_stats_int_covariance_m(c2.data(), 1, d2.data(), 1, c2.size(), 3, 4);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::cov(c + c2.cast<double>(), d + d2.cast<double>());
}

TEST_CASE("stat_corrcoef")
{
    double v, g;
    double work[20];

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::corrcoef(c, d);
    g = gsl_stats_correlation(c.data(), 1, d.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2, d2;
    c2 = iexp::ArrayXi::Random(10);
    d2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::corrcoef(c2, d2);
    g = gsl_stats_int_correlation(c2.data(), 1, d2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::spearman(c, d);
    g = gsl_stats_spearman(c.data(), 1, d.data(), 1, c.size(), work);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::ArrayXi::Random(10);
    d2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::spearman(c2, d2);
    g = gsl_stats_int_spearman(c2.data(), 1, d2.data(), 1, c2.size(), work);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::corrcoef(c + c2.cast<double>(), d + d2.cast<double>());
    v = iexp::stats::spearman(c + c2.cast<double>(), d + d2.cast<double>());
}
