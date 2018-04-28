#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <sort/sort.h>
#include <stats/autocorr.h>
#include <stats/corrcoef.h>
#include <stats/cov.h>
#include <stats/kurtosis.h>
#include <stats/mean.h>
#include <stats/median.h>
#include <stats/quantile.h>
#include <stats/skewness.h>
#include <stats/std.h>
#include <stats/tss.h>
#include <stats/var.h>

TEST_CASE("stat_mean")
{
    double v, g;

    iexp::ArrayXd c, w;
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

    c = iexp::ArrayXd::Random(10);
    w = iexp::ArrayXd::Random(10);
    v = iexp::stats::wmean(c, w);
    g = gsl_stats_wmean(w.data(), 1, c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::wmean(c + c2.cast<double>(), c + c2.cast<double>());
}

TEST_CASE("stat_var")
{
    double v, g;

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::var(c);
    g = gsl_stats_variance(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wvar(c, d);
    g = gsl_stats_wvariance(d.data(), 1, c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::var(c, 1);
    g = gsl_stats_variance_m(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wvar(c, d, 2);
    g = gsl_stats_wvariance_m(d.data(), 1, c.data(), 1, c.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::uvar(c, 1);
    g = gsl_stats_variance_with_fixed_mean(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wuvar(c, d, 2);
    g = gsl_stats_wvariance_with_fixed_mean(d.data(),
                                            1,
                                            c.data(),
                                            1,
                                            c.size(),
                                            2);
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

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::std(c);
    g = gsl_stats_sd(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wstd(c, d);
    g = gsl_stats_wsd(d.data(), 1, c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::std(c, 1);
    g = gsl_stats_sd_m(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wstd(c, d, 2);
    g = gsl_stats_wsd_m(d.data(), 1, c.data(), 1, c.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::ustd(c, 1);
    g = gsl_stats_sd_with_fixed_mean(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wustd(c, d, 2);
    g = gsl_stats_wsd_with_fixed_mean(d.data(), 1, c.data(), 1, c.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::abstd(c);
    g = gsl_stats_absdev(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wabstd(c, d);
    g = gsl_stats_wabsdev(d.data(), 1, c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::abstd(c, 2);
    g = gsl_stats_absdev_m(c.data(), 1, c.size(), 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wabstd(c, d, 2);
    g = gsl_stats_wabsdev_m(d.data(), 1, c.data(), 1, c.size(), 2);
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

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::tss(c);
    g = gsl_stats_tss(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wtss(c, d);
    g = gsl_stats_wtss(d.data(), 1, c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::tss(c, 1);
    g = gsl_stats_tss_m(c.data(), 1, c.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wtss(c, d, 2);
    g = gsl_stats_wtss_m(d.data(), 1, c.data(), 1, c.size(), 2);
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

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::skewness(c);
    g = gsl_stats_skew(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wskewness(c, d);
    g = gsl_stats_wskew(d.data(), 1, c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::skewness(c, 1, 2);
    g = gsl_stats_skew_m_sd(c.data(), 1, c.size(), 1, 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wskewness(c, d, 2, 3);
    g = gsl_stats_wskew_m_sd(d.data(), 1, c.data(), 1, c.size(), 2, 3);
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

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::kurtosis(c);
    g = gsl_stats_kurtosis(c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wkurtosis(c, d);
    g = gsl_stats_wkurtosis(d.data(), 1, c.data(), 1, c.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::kurtosis(c, 1, 2);
    g = gsl_stats_kurtosis_m_sd(c.data(), 1, c.size(), 1, 2);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c = iexp::ArrayXd::Random(10);
    d = iexp::ArrayXd::Random(10);
    v = iexp::stats::wkurtosis(c, d, 3, 4);
    g = gsl_stats_wkurtosis_m_sd(d.data(), 1, c.data(), 1, c.size(), 3, 4);
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

    iexp::MatrixXi c2(3, 4);
    c2 = iexp::MatrixXi::Random(3, 4);
    v = iexp::stats::autocorr(c2);
    g = gsl_stats_int_lag1_autocorrelation(c2.data(), 1, c2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    c2 = iexp::MatrixXi::Random(3, 4);
    v = iexp::stats::autocorr(c2, 1);
    g = gsl_stats_int_lag1_autocorrelation_m(c2.data(), 1, c2.size(), 1);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    // compile
    v = iexp::stats::mean(c2.cast<double>() + c2.cast<double>());
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

TEST_CASE("stat_median")
{
    double v, g;
    double work[20];

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::median(iexp::sort(c));
    d = sort(c);
    g = gsl_stats_median_from_sorted_data(d.data(), 1, d.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2, d2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::median(iexp::sort(c2));
    d2 = sort(c2);
    g = gsl_stats_int_median_from_sorted_data(d2.data(), 1, d2.size());
    REQUIRE(__D_EQ_IN(v, g, 1e-9));
}

TEST_CASE("stat_qunatile")
{
    double v, g;
    double work[20];

    iexp::ArrayXd c, d;
    c = iexp::ArrayXd::Random(10);
    v = iexp::stats::quantile(iexp::sort(c), 0.93);
    d = sort(c);
    g = gsl_stats_quantile_from_sorted_data(d.data(), 1, d.size(), 0.93);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));

    iexp::ArrayXi c2, d2;
    c2 = iexp::ArrayXi::Random(10);
    v = iexp::stats::quantile(iexp::sort(c2), 0.88);
    d2 = sort(c2);
    g = gsl_stats_int_quantile_from_sorted_data(d2.data(), 1, d2.size(), 0.88);
    REQUIRE(__D_EQ_IN(v, g, 1e-9));
}
