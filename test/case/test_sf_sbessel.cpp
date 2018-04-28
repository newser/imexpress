#include <catch.hpp>
#include <iostream>
#include <special/modified_spherical_bessel.h>
#include <special/spherical_bessel.h>
#include <test_util.h>

TEST_CASE("sf_sbessel_j")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // J0
    m << -10.0, 0.001, 1.0, 100.0;
    m2 = iexp::sf::sbessel_j<iexp::sf::j0>(m);
    REQUIRE(__D_EQ9(m2(0, 0), -0.05440211108893698134));
    REQUIRE(__D_EQ9(m2(1, 0), 0.9999998333333416667));
    REQUIRE(__D_EQ9(m2(2, 0), 0.84147098480789650670));
    REQUIRE(__D_EQ9(m2(3, 0), -0.005063656411097587937));

    m2 = iexp::sf::sbessel_j<iexp::sf::j0>(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -0.05440211108893698134));
    REQUIRE(__D_EQ9(m2(1, 0), 0.9999998333333416667));
    REQUIRE(__D_EQ9(m2(2, 0), 0.84147098480789650670));
    REQUIRE(__D_EQ9(m2(3, 0), -0.005063656411097587937));
    gsl_sf_bessel_j0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_j0_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    __draw_1d(sbessel_j<iexp::sf::j0>, 0, 10, 100, 0, 0);

    // J1
    m << -10.0, 0.01, 1.0, 100.0;
    m2 = iexp::sf::sbessel_j<iexp::sf::j1>(m);
    REQUIRE(__D_EQ9(m2(0, 0), -0.07846694179875154709));
    REQUIRE(__D_EQ9(m2(1, 0), 0.003333300000119047399));
    REQUIRE(__D_EQ9(m2(2, 0), 0.30116867893975678925));
    REQUIRE(__D_EQ9(m2(3, 0), -0.008673825286987815220));

    m2 = iexp::sf::sbessel_j<iexp::sf::j1>(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -0.07846694179875154709));
    REQUIRE(__D_EQ9(m2(1, 0), 0.003333300000119047399));
    REQUIRE(__D_EQ9(m2(2, 0), 0.30116867893975678925));
    REQUIRE(__D_EQ9(m2(3, 0), -0.008673825286987815220));
    gsl_sf_bessel_j1_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_j1_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    __draw_1d(sbessel_j<iexp::sf::j1>, 0, 10, 100, 0, 0);

    // J2
    m2 = iexp::sf::sbessel_j<iexp::sf::j2>(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.07794219362856244547));
    REQUIRE(__D_EQ9(m2(1, 0), 6.666619047751322551e-06));
    REQUIRE(__D_EQ9(m2(2, 0), 0.06203505201137386110));
    REQUIRE(__D_EQ9(m2(3, 0), 0.004803441652487953480));

    m2 = iexp::sf::sbessel_j<iexp::sf::j2>(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.07794219362856244547));
    REQUIRE(__D_EQ9(m2(1, 0), 6.666619047751322551e-06));
    REQUIRE(__D_EQ9(m2(2, 0), 0.06203505201137386110));
    REQUIRE(__D_EQ9(m2(3, 0), 0.004803441652487953480));
    gsl_sf_bessel_j2_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_j2_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    __draw_1d(sbessel_j<iexp::sf::j2>, 0, 10, 100, 0, 0);

    // Jn
    m << 1.0, 5.001, 10.0, 100.0;
    m2 = iexp::sf::sbessel_j(100, m);
    REQUIRE(__D_EQ9(m2(3, 0), 0.010880477011438336539));

    m2 = iexp::sf::sbessel_j(100, m, e);
    REQUIRE(__D_EQ9(m2(3, 0), 0.010880477011438336539));
    gsl_sf_bessel_jl_e(100, m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));
}

TEST_CASE("sf_sbessel_y")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // y0
    m << 0.001, 1.0, 100.0, 65536.0;
    m2 = iexp::sf::sbessel_y<iexp::sf::y0>(m);
    REQUIRE(__D_EQ9(m2(0, 0), -999.99950000004166670));
    REQUIRE(__D_EQ9(m2(1, 0), -0.5403023058681397174));
    REQUIRE(__D_EQ9(m2(2, 0), -0.008623188722876839341));
    REQUIRE(__D_EQ9(m2(3, 0), 0.000011014324202158573930));

    m2 = iexp::sf::sbessel_y<iexp::sf::y0>(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -999.99950000004166670));
    REQUIRE(__D_EQ9(m2(3, 0), 0.000011014324202158573930));
    gsl_sf_bessel_y0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_y0_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    __draw_1d(sbessel_y<iexp::sf::y0>, 0.01, 0.1, 100, 0, 0);

    // y1
    m << 0.01, 1.0, 10.0, 100.0;
    m2 = iexp::sf::sbessel_y<iexp::sf::y1>(m);
    REQUIRE(__D_EQ9(m2(0, 0), -10000.499987500069444));
    REQUIRE(__D_EQ9(m2(1, 0), -1.3817732906760362241));
    REQUIRE(__D_EQ9(m2(2, 0), 0.06279282637970150586));
    REQUIRE(__D_EQ9(m2(3, 0), 0.004977424523868819543));

    m2 = iexp::sf::sbessel_y<iexp::sf::y1>(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -10000.499987500069444));
    REQUIRE(__D_EQ9(m2(3, 0), 0.004977424523868819543));
    gsl_sf_bessel_y1_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_y1_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    __draw_1d(sbessel_y<iexp::sf::y1>, 0.01, 0.1, 100, 0, 0);

    // y2
    m << 0.01, 1.0, 10.0, 100.0;
    m2 = iexp::sf::sbessel_y<iexp::sf::y2>(m);
    REQUIRE(__D_EQ9(m2(0, 0), -3.0000500012499791668e+06));
    REQUIRE(__D_EQ9(m2(1, 0), -3.605017566159968955));
    REQUIRE(__D_EQ9(m2(2, 0), -0.06506930499373479347));
    REQUIRE(__D_EQ9(m2(3, 0), 0.008772511458592903927));

    m2 = iexp::sf::sbessel_y<iexp::sf::y2>(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -3.0000500012499791668e+06));
    REQUIRE(__D_EQ9(m2(3, 0), 0.008772511458592903927));
    gsl_sf_bessel_y2_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_y2_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    __draw_1d(sbessel_y<iexp::sf::y2>, 0.01, 0.1, 100, 0, 0);

    // yl
    m.resize(2);
    m << 0.01, 10.0;
    m2 = iexp::sf::sbessel_y(10, m);
    // REQUIRE(__D_EQ9(m2(0, 0), -6.5473079797378378e+30));
    REQUIRE(__D_EQ9(m2(1, 0), -0.172453672088057849));

    m2 = iexp::sf::sbessel_y(10, m, e);
    // REQUIRE(__D_EQ9(m2(0, 0), -6.5473079797378378e+30));
    REQUIRE(__D_EQ9(m2(1, 0), -0.172453672088057849));
    gsl_sf_bessel_yl_e(10, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_yl_e(10, m(1, 0), &r);
    REQUIRE(__D_EQ9(e(1, 0), r.err));
}

TEST_CASE("sf_sbessel_i_scaled")
{
    iexp::ArrayXd m(3), m2(3), e(3);
    gsl_sf_result r;

    // I0
    m << 0.1, 2.0, 100.0;
    m2 = iexp::sf::sbessel_i<true>(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.9063462346100907067));
    REQUIRE(__D_EQ9(m2(1, 0), 0.24542109027781645493));
    REQUIRE(__D_EQ9(m2(2, 0), 0.005000000000000000000));

    m2 = iexp::sf::sbessel_i<true>(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.9063462346100907067));
    REQUIRE(__D_EQ9(m2(1, 0), 0.24542109027781645493));
    REQUIRE(__D_EQ9(m2(2, 0), 0.005000000000000000000));
    gsl_sf_bessel_i0_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_i0_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ7(e(2, 0), r.err));

    // I1
    m2 = iexp::sf::sbessel_i<true>(1, m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.030191419289002226846));
    REQUIRE(__D_EQ9(m2(1, 0), 0.131868364583275317610));
    REQUIRE(__D_EQ9(m2(2, 0), 0.004950000000000000000));

    m2 = iexp::sf::sbessel_i<true>(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.030191419289002226846));
    REQUIRE(__D_EQ9(m2(1, 0), 0.131868364583275317610));
    REQUIRE(__D_EQ9(m2(2, 0), 0.004950000000000000000));
    gsl_sf_bessel_i1_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_i1_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // In
    m2 = iexp::sf::sbessel_i<true>(4, m);
    REQUIRE(__D_EQ9(m2(0, 0), 9.579352242057134927e-08));

    m2 = iexp::sf::sbessel_i<true>(4, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 9.579352242057134927e-08));
    gsl_sf_bessel_In_e(4, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
}

TEST_CASE("sf_sbessel_k")
{
    iexp::ArrayXd m(3), m2(3), e(3);
    gsl_sf_result r;

    // K0
    m << 0.1, 2.0, 100.0;
    m2 = iexp::sf::sbessel_k<true>(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), 15.707963267948966192));
    REQUIRE(__D_EQ9(m2(1, 0), 0.7853981633974483096));
    REQUIRE(__D_EQ9(m2(2, 0), 0.015707963267948966192));

    m2 = iexp::sf::sbessel_k<true>(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 15.707963267948966192));
    REQUIRE(__D_EQ9(m2(1, 0), 0.7853981633974483096));
    REQUIRE(__D_EQ9(m2(2, 0), 0.015707963267948966192));
    gsl_sf_bessel_k0_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_k0_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // K1
    m << 0.1, 2.0, 100.0;
    m2 = iexp::sf::sbessel_k<true>(1, m);
    // std::cout << m2 << std::endl;
    REQUIRE(__D_EQ9(m2(0, 0), 172.78759594743862812));
    REQUIRE(__D_EQ9(m2(1, 0), 1.1780972450961724644));
    REQUIRE(__D_EQ9(m2(2, 0), 0.015865042900628455854));

    m2 = iexp::sf::sbessel_k<true>(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 172.78759594743862812));
    REQUIRE(__D_EQ9(m2(1, 0), 1.1780972450961724644));
    REQUIRE(__D_EQ9(m2(2, 0), 0.015865042900628455854));
    gsl_sf_bessel_k1_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_k1_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // Kn
    m.resize(2);
    m << 1.0 / 256.0, 1.0 / 8.0;
    m2 = iexp::sf::sbessel_k<true>(4, m);
    // REQUIRE(__D_EQ_Nep(m2(0, 0), 1.8205599816961954439e+14, 1000));
    REQUIRE(__D_EQ9(m2(1, 0), 6.1173217814406597530e+06));

    m2 = iexp::sf::sbessel_k<true>(4, m, e);
    // REQUIRE(__D_EQ9(m2(0, 0), 1.8205599816961954439e+14));
    REQUIRE(__D_EQ9(m2(1, 0), 6.1173217814406597530e+06));
    gsl_sf_bessel_kl_scaled_e(4, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_kl_scaled_e(4, m(1, 0), &r);
    REQUIRE(__D_EQ9(e(1, 0), r.err));
}
