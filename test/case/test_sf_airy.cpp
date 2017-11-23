#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <special/airy.h>
#include <special/airy_zero.h>

TEST_CASE("sf_airy")
{
    iexp::ArrayXXd x = iexp::ArrayXXd::Random(2, 2), y(2, 2), e(2, 2);

    x << -500.0, -5.0, -0.3000000000000094, 0.6999999999999907;
    y = iexp::sf::airy_Ai(x);
    REQUIRE(__D_EQ9(y(0, 0), 0.0725901201040411396));
    REQUIRE(__D_EQ9(y(0, 1), 0.3507610090241142));
    REQUIRE(__D_EQ9(y(1, 0), 0.4309030952855831));
    REQUIRE(__D_EQ9(y(1, 1), 0.1891624003981519));

    y = iexp::sf::airy_Bi(x);
    REQUIRE(__D_EQ9(y(0, 0), -0.094688570132991028));
    REQUIRE(__D_EQ9(y(0, 1), -0.1383691349016005));
    REQUIRE(__D_EQ9(y(1, 1), 0.9733286558781599));

    y = iexp::sf::airy_Ai(x, e);
    REQUIRE(__D_EQ9(y(0, 0), 0.0725901201040411396));
    REQUIRE(__D_EQ9(y(0, 1), 0.3507610090241142));
    REQUIRE(__D_EQ9(y(1, 0), 0.4309030952855831));
    REQUIRE(__D_EQ9(y(1, 1), 0.1891624003981519));

    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_e(-500.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Ai_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_Ai_e(-0.3000000000000094, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
        gsl_sf_airy_Ai_e(0.6999999999999907, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 1), r.err));
    }

    y = iexp::sf::airy_Bi(x, e);
    REQUIRE(__D_EQ9(y(0, 0), -0.094688570132991028));
    REQUIRE(__D_EQ9(y(0, 1), -0.1383691349016005));
    REQUIRE(__D_EQ9(y(1, 1), 0.9733286558781599));

    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_e(-500.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Bi_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_Bi_e(-0.3000000000000094, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
        gsl_sf_airy_Bi_e(0.6999999999999907, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 1), r.err));
    }

    // compile
    y = iexp::sf::airy_Ai((x + x) * 2);
    y = iexp::sf::airy_Ai((x + x) * 2, e);
    // y = iexp::sf::airy_Ai((x + x) * 2, e + e);
    y = iexp::sf::airy_Bi((x + x) * 2);
    y = iexp::sf::airy_Bi((x + x) * 2, e);
}

TEST_CASE("sf_airy_scaled")
{
    iexp::ArrayXXd x = iexp::ArrayXXd::Random(2, 3), y(2, 3), e(2, 3);

    x << -5.0, 0.6999999999999907, 1.649999999999991, 2.54999999999999,
        3.499999999999987, 5.39999999999998;
    y = iexp::sf::airy_Ai<true>(x);
    REQUIRE(__D_EQ9(y(0, 0), 0.3507610090241142));
    REQUIRE(__D_EQ9(y(0, 1), 0.2795125667681217));
    REQUIRE(__D_EQ9(y(0, 2), 0.2395493001442741));
    REQUIRE(__D_EQ9(y(1, 0), 0.2183658595899388));
    REQUIRE(__D_EQ9(y(1, 1), 0.2032920808163519));
    REQUIRE(__D_EQ9(y(1, 2), 0.1836050093282229));

    y = iexp::sf::airy_Bi<true>(x);
    REQUIRE(__D_EQ9(y(0, 0), -0.1383691349016005));
    REQUIRE(__D_EQ9(y(1, 2), 0.3734050675720473));

    y = iexp::sf::airy_Ai<true>(x, e);
    REQUIRE(__D_EQ9(y(0, 0), 0.3507610090241142));
    REQUIRE(__D_EQ9(y(0, 1), 0.2795125667681217));
    REQUIRE(__D_EQ9(y(0, 2), 0.2395493001442741));
    REQUIRE(__D_EQ9(y(1, 0), 0.2183658595899388));
    REQUIRE(__D_EQ9(y(1, 1), 0.2032920808163519));
    REQUIRE(__D_EQ9(y(1, 2), 0.1836050093282229));

    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_scaled_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Ai_scaled_e(0.6999999999999907, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_Ai_scaled_e(5.39999999999998, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
    }

    y = iexp::sf::airy_Bi<true>(x, e);
    REQUIRE(__D_EQ9(y(0, 0), -0.1383691349016005));
    REQUIRE(__D_EQ9(y(1, 2), 0.3734050675720473));

    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_scaled_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Bi_scaled_e(0.6999999999999907, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_Bi_scaled_e(5.39999999999998, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
    }

    // compile
    y = iexp::sf::airy_Ai<true>((x + x) * 2);
    y = iexp::sf::airy_Ai<true>((x + x) * 2, e);
    // y = iexp::sf::airy_Ai((x + x) * 2, e + e);
    y = iexp::sf::airy_Bi<true>((x + x) * 2);
    y = iexp::sf::airy_Bi<true>((x + x) * 2, e);
}

TEST_CASE("sf_airy_deriv")
{
    iexp::ArrayXXd x = iexp::ArrayXXd::Random(2, 2), y(2, 2), e(2, 2);

    x << -5.0, -0.5500000000000094, 0.4999999999999906, 5.199999999999981;
    y = iexp::sf::airy_Ai_deriv(x);
    REQUIRE(__D_EQ9(y(0, 0), 0.3271928185544435));
    REQUIRE(__D_EQ9(y(0, 1), -0.1914604987143629));
    REQUIRE(__D_EQ9(y(1, 0), -0.2249105326646850));
    REQUIRE(__D_EQ9(y(1, 1), -0.0001589434526459543));

    y = iexp::sf::airy_Bi_deriv(x);
    REQUIRE(__D_EQ9(y(0, 0), 0.778411773001899));
    REQUIRE(__D_EQ9(y(0, 1), 0.5155785358765014));
    REQUIRE(__D_EQ9(y(1, 0), 0.5445725641405883));
    REQUIRE(__D_EQ9(y(1, 1), 2279.748293583233));

    y = iexp::sf::airy_Ai_deriv(x, e);
    REQUIRE(__D_EQ9(y(0, 0), 0.3271928185544435));
    REQUIRE(__D_EQ9(y(0, 1), -0.1914604987143629));
    REQUIRE(__D_EQ9(y(1, 0), -0.2249105326646850));
    REQUIRE(__D_EQ9(y(1, 1), -0.0001589434526459543));

    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Ai_e(-0.5500000000000094, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_Ai_e(0.4999999999999906, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
        gsl_sf_airy_Ai_e(5.199999999999981, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 1), r.err));
    }

    y = iexp::sf::airy_Bi_deriv(x, e);
    REQUIRE(__D_EQ9(y(0, 0), 0.778411773001899));
    REQUIRE(__D_EQ9(y(0, 1), 0.5155785358765014));
    REQUIRE(__D_EQ9(y(1, 0), 0.5445725641405883));
    REQUIRE(__D_EQ9(y(1, 1), 2279.748293583233));

    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Bi_e(-0.5500000000000094, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_Bi_e(0.4999999999999906, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
        gsl_sf_airy_Bi_e(5.199999999999981, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 1), r.err));
    }

    // compile
    y = iexp::sf::airy_Ai((x + x) * 2);
    y = iexp::sf::airy_Ai((x + x) * 2, e);
    // y = iexp::sf::airy_Ai((x + x) * 2, e + e);
    y = iexp::sf::airy_Bi((x + x) * 2);
    y = iexp::sf::airy_Bi((x + x) * 2, e);
}

TEST_CASE("sf_airy_deriv_scaled")
{
    iexp::ArrayXXd x = iexp::ArrayXXd::Random(2, 1), y(2, 1), e(2, 1);

    x << -5.0, 6.299999999999977;
    y = iexp::sf::airy_Ai_deriv<true>(x);
    REQUIRE(__D_EQ9(y(0, 0), 0.3271928185544435));
    REQUIRE(__D_EQ9(y(1, 0), -0.4508799189585947));

    y = iexp::sf::airy_Bi_deriv<true>(x);
    REQUIRE(__D_EQ9(y(0, 0), 0.778411773001899));
    REQUIRE(__D_EQ9(y(1, 0), 0.8852064139737571));

    y = iexp::sf::airy_Ai_deriv<true>(x, e);
    REQUIRE(__D_EQ9(y(0, 0), 0.3271928185544435));
    REQUIRE(__D_EQ9(y(1, 0), -0.4508799189585947));

    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_deriv_scaled_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Ai_deriv_scaled_e(6.299999999999977, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
    }

    y = iexp::sf::airy_Bi_deriv<true>(x, e);
    REQUIRE(__D_EQ9(y(0, 0), 0.778411773001899));
    REQUIRE(__D_EQ9(y(1, 0), 0.8852064139737571));

    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_deriv_scaled_e(-5.0, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_Bi_deriv_scaled_e(6.299999999999977, GSL_PREC_DOUBLE, &r);
        REQUIRE(__D_EQ9(e(1, 0), r.err));
    }

    // compile
    y = iexp::sf::airy_Ai_deriv<true>((x + x) * 2);
    y = iexp::sf::airy_Ai_deriv<true>((x + x) * 2, e);
    // y = iexp::sf::airy_Ai((x + x) * 2, e + e);
}

TEST_CASE("sf_airy_zero")
{
    iexp::ArrayXXi x = iexp::ArrayXXi::Random(1, 4);
    iexp::ArrayXXd y(1, 4), e(1, 4);

    x << 2, 50, 100, 110;
    y = iexp::sf::airy_n0_Ai(x);
    REQUIRE(__D_EQ9(y(0, 0), -4.087949444130970617));
    REQUIRE(__D_EQ9(y(0, 1), -38.02100867725525443));
    REQUIRE(__D_EQ9(y(0, 2), -60.45555727411669871));
    REQUIRE(__D_EQ9(y(0, 3), -64.43135670991324811));

    y = iexp::sf::airy_n0_Bi(x);
    REQUIRE(__D_EQ9(y(0, 0), -3.271093302836352716));
    REQUIRE(__D_EQ9(y(0, 1), -37.76583438165180116));
    REQUIRE(__D_EQ9(y(0, 2), -60.25336482580837088));
    REQUIRE(__D_EQ9(y(0, 3), -64.2355167606561537));

    y = iexp::sf::airy_n0_Ai(x, e);
    REQUIRE(__D_EQ9(y(0, 0), -4.087949444130970617));
    REQUIRE(__D_EQ9(y(0, 1), -38.02100867725525443));
    REQUIRE(__D_EQ9(y(0, 2), -60.45555727411669871));
    REQUIRE(__D_EQ9(y(0, 3), -64.43135670991324811));

    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Ai(2);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_zero_Ai(50);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_zero_Ai(100);
        REQUIRE(__D_EQ9(e(0, 2), r.err));
        gsl_sf_airy_zero_Ai(110);
        REQUIRE(__D_EQ9(e(0, 3), r.err));
    }

    y = iexp::sf::airy_n0_Bi(x, e);
    REQUIRE(__D_EQ9(y(0, 0), -3.271093302836352716));
    REQUIRE(__D_EQ9(y(0, 1), -37.76583438165180116));
    REQUIRE(__D_EQ9(y(0, 2), -60.25336482580837088));
    REQUIRE(__D_EQ9(y(0, 3), -64.2355167606561537));

    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Bi(2);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_zero_Bi(50);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_zero_Bi(100);
        REQUIRE(__D_EQ9(e(0, 2), r.err));
        gsl_sf_airy_zero_Bi(110);
        REQUIRE(__D_EQ9(e(0, 3), r.err));
    }
}

TEST_CASE("sf_airy_zero_deriv")
{
    iexp::ArrayXXi x = iexp::ArrayXXi::Random(1, 4);
    iexp::ArrayXXd y(1, 4), e(1, 4);

    x << 2, 50, 100, 110;
    y = iexp::sf::airy_n0_Ai_deriv(x);
    REQUIRE(__D_EQ9(y(0, 0), -3.248197582179836561));
    REQUIRE(__D_EQ9(y(0, 1), -37.76565910053887108));
    REQUIRE(__D_EQ9(y(0, 2), -60.25329596442479317));
    REQUIRE(__D_EQ9(y(0, 3), -64.23545617243546956));

    y = iexp::sf::airy_n0_Bi_deriv(x);
    REQUIRE(__D_EQ9(y(0, 0), -4.073155089071828216));
    REQUIRE(__D_EQ9(y(0, 1), -38.02083574095788210));
    REQUIRE(__D_EQ9(y(0, 2), -60.45548887257140819));
    REQUIRE(__D_EQ9(y(0, 3), -64.43129648944845060));

    y = iexp::sf::airy_n0_Ai_deriv(x, e);
    REQUIRE(__D_EQ9(y(0, 0), -3.248197582179836561));
    REQUIRE(__D_EQ9(y(0, 1), -37.76565910053887108));
    REQUIRE(__D_EQ9(y(0, 2), -60.25329596442479317));
    REQUIRE(__D_EQ9(y(0, 3), -64.23545617243546956));

    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Ai_deriv(2);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_zero_Ai_deriv(50);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_zero_Ai_deriv(100);
        REQUIRE(__D_EQ9(e(0, 2), r.err));
        gsl_sf_airy_zero_Ai_deriv(110);
        REQUIRE(__D_EQ9(e(0, 3), r.err));
    }

    y = iexp::sf::airy_n0_Bi_deriv(x, e);
    REQUIRE(__D_EQ9(y(0, 0), -4.073155089071828216));
    REQUIRE(__D_EQ9(y(0, 1), -38.02083574095788210));
    REQUIRE(__D_EQ9(y(0, 2), -60.45548887257140819));
    REQUIRE(__D_EQ9(y(0, 3), -64.43129648944845060));

    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Bi_deriv(2);
        REQUIRE(__D_EQ9(e(0, 0), r.err));
        gsl_sf_airy_zero_Bi_deriv(50);
        REQUIRE(__D_EQ9(e(0, 1), r.err));
        gsl_sf_airy_zero_Bi_deriv(100);
        REQUIRE(__D_EQ9(e(0, 2), r.err));
        gsl_sf_airy_zero_Bi_deriv(110);
        REQUIRE(__D_EQ9(e(0, 3), r.err));
    }
}
