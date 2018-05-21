#include <catch.hpp>
#include <iostream>
#include <special/psi.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_psi")
{
    Matrix2i p;
    Matrix2d r, e;
    p << 1, 2, 3, 4;
    r = sf::psi(p);
    REQUIRE(__D_EQ9(r(0, 0), -0.57721566490153286060));
    REQUIRE(__D_EQ9(r(0, 1), 0.42278433509846713939));
    REQUIRE(__D_EQ9(r(1, 0), 0.92278433509846713939));
    REQUIRE(__D_EQ9(r(1, 1), 1.2561176684318004727));

    __draw_1d(psi, 0.1, 10, 100, 0, 0);

    r = sf::psi(p, e);
    REQUIRE(__D_EQ9(r(0, 0), -0.57721566490153286060));
    REQUIRE(__D_EQ9(r(0, 1), 0.42278433509846713939));
    REQUIRE(__D_EQ9(r(1, 0), 0.92278433509846713939));
    REQUIRE(__D_EQ9(r(1, 1), 1.2561176684318004727));

    gsl_sf_result re;
    gsl_sf_psi_int_e(1, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_psi_int_e(4, &re);
    REQUIRE(e(1, 1) == re.err);

    // double
    Matrix2d p2;
    p2 << -100.5, -10.5, 5.0, 5000.0;
    r = sf::psi(p2);
    REQUIRE(__D_EQ9(r(0, 0), 4.615124601338064117));
    REQUIRE(__D_EQ9(r(0, 1), 2.3982391295357816134));
    REQUIRE(__D_EQ9(r(1, 0), 1.5061176684318004727));
    REQUIRE(__D_EQ9(r(1, 1), 8.517093188082904107));

    r = sf::psi(p2, e);
    REQUIRE(__D_EQ9(r(0, 0), 4.615124601338064117));
    REQUIRE(__D_EQ9(r(0, 1), 2.3982391295357816134));
    REQUIRE(__D_EQ9(r(1, 0), 1.5061176684318004727));
    REQUIRE(__D_EQ9(r(1, 1), 8.517093188082904107));

    gsl_sf_psi_e(-100.5, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_psi_e(5000.0, &re);
    REQUIRE(e(1, 1) == re.err);

    // 1pix
    p2 << 0.8, 1.0, 5.0, 100.0;
    r = sf::psi_1pix(p2);
    REQUIRE(__D_EQ9(r(0, 0), -0.07088340212750589223));
    REQUIRE(__D_EQ9(r(1, 1), 4.605178519404762003));

    r = sf::psi_1pix(p2, e);
    REQUIRE(__D_EQ9(r(0, 0), -0.07088340212750589223));
    REQUIRE(__D_EQ9(r(1, 1), 4.605178519404762003));
    gsl_sf_psi_1piy_e(0.8, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_psi_1piy_e(100.0, &re);
    REQUIRE(e(1, 1) == re.err);
}

TEST_CASE("test_psi_d1")
{
    Matrix2i p;
    Matrix2d r, e;
    p << 1, 2, 3, 4;
    r = sf::psi_d1(p);
    REQUIRE(__D_EQ9(r(0, 0), 1.6449340668482264364));
    REQUIRE(__D_EQ9(r(0, 1), 0.64493406684822643647));
    REQUIRE(__D_EQ9(r(1, 0), 0.39493406684822643647));
    REQUIRE(__D_EQ9(r(1, 1), 0.28382295573711532536));

    __draw_1d(psi_d1, -10, 10, 100, 0, 0);

    r = sf::psi_d1(p, e);
    REQUIRE(__D_EQ9(r(0, 0), 1.6449340668482264364));
    REQUIRE(__D_EQ9(r(0, 1), 0.64493406684822643647));
    REQUIRE(__D_EQ9(r(1, 0), 0.39493406684822643647));
    REQUIRE(__D_EQ9(r(1, 1), 0.28382295573711532536));

    gsl_sf_result re;
    gsl_sf_psi_1_int_e(1, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_psi_1_int_e(4, &re);
    REQUIRE(e(1, 1) == re.err);

    // double
    Matrix2d p2;
    p2 << -1.0 - 1.0 / 128.0, -1.50, 1.0 / 32.0, 500.0;
    r = sf::psi_d1(p2);
    REQUIRE(__D_EQ9(r(0, 0), 16386.648472598746587));
    REQUIRE(__D_EQ9(r(0, 1), 9.3792466449891237539));
    REQUIRE(__D_EQ9(r(1, 0), 1025.5728544782377089));
    REQUIRE(__D_EQ9(r(1, 1), 0.0020020013333322666697));

    r = sf::psi_d1(p2, e);
    REQUIRE(__D_EQ9(r(0, 0), 16386.648472598746587));
    REQUIRE(__D_EQ9(r(0, 1), 9.3792466449891237539));
    REQUIRE(__D_EQ9(r(1, 0), 1025.5728544782377089));
    REQUIRE(__D_EQ9(r(1, 1), 0.0020020013333322666697));

    gsl_sf_psi_1_e(-1.0 - 1.0 / 128.0, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_psi_1_e(500.0, &re);
    REQUIRE(e(1, 1) == re.err);
}

TEST_CASE("test_psi_dn")
{
    Vector2d p, r, e;
    p << 5.0, 500.0;
    r = sf::psi_dn(3, p);
    REQUIRE(__D_EQ9(r(0), 0.021427828192755075022));
    REQUIRE(__D_EQ9(r(1), 1.6048063999872000683e-08));

    r = sf::psi_dn(3, p, e);
    REQUIRE(__D_EQ9(r(0), 0.021427828192755075022));
    REQUIRE(__D_EQ9(r(1), 1.6048063999872000683e-08));

    gsl_sf_result re;
    gsl_sf_psi_n_e(3, 5.0, &re);
    REQUIRE(e(0) == re.err);
    gsl_sf_psi_n_e(3, 500.0, &re);
    REQUIRE(e(1) == re.err);
}
